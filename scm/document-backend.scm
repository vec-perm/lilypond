;;; backend-documentation-lib.scm -- Functions for backend documentation
;;;
;;; source file of the GNU LilyPond music typesetter
;;; 
;;; (c)  2000--2004 Han-Wen Nienhuys <hanwen@cs.uu.nl>
;;; Jan Nieuwenhuizen <janneke@gnu.org>


;;; This file generates documentation for the backend of lilypond.

;; alist of property descriptions

;;
"
TODO:


Grob bla

Created by:

  * preset properties + explanation

Interfaces:

  * properties available.

"


(define (interface-doc-string interface grob-description)
  (let*
      ((name (car interface))
       (desc (cadr interface))
       (props (sort (caddr interface) symbol<?))
       (docfunc (lambda (pr)
		  (property->texi
		    'backend pr grob-description)))
       (iprops (filter (lambda (x) (object-property x 'backend-internal) ) props))
       (uprops (filter (lambda (x) (not (object-property x 'backend-internal)) ) props))
       (user-propdocs (map docfunc uprops))
       (internal-propdocs (map docfunc iprops)))

       (string-append
	desc
	"\n\n@unnumberedsubsubsec User settable properties:\n"
	(description-list->texi user-propdocs)

	"\n\n@unnumberedsubsubsec Internal properties: \n"
	(description-list->texi internal-propdocs)
	)
    ))


(define iface->grob-table (make-vector 61 '()))
;; extract ifaces, and put grob into the hash table.
(map
 (lambda (x)
   (let*
       (
	(metah (assoc 'meta (cdr x)))
	(meta (cdr metah))
	(ifaces (cdr (assoc 'interfaces meta)))
	)

     (map (lambda (iface)
	    (hashq-set!
	     iface->grob-table iface
	     (cons (car x)
		   (hashq-ref iface->grob-table iface '())
		   )))
	  ifaces)
     ))
 all-grob-descriptions)

;; First level Interface description
(define (interface-doc interface)
  (let ((name (symbol->string (car interface))))
    (make <texi-node>
      #:name name
      #:text (string-append
	      (interface-doc-string (cdr interface) '())
	      "\n\n"
	      "This grob interface is used in the following graphical objects: "

	      (human-listify
	       (map ref-ify
		    (map symbol->string
			 (hashq-ref iface->grob-table (car interface) '() )))))

      )))

(define (grob-alist->texi alist)
  (let*
      ((uprops (filter (lambda (x) (not (object-property x 'backend-internal)))
		       (map car alist))))

    (description-list->texi
     (map (lambda (y) (property->texi 'backend y alist))
	  uprops)
     )))


(define (grob-doc description)
  "Given a property alist DESCRIPTION, make a documentation
node."
  
  (let*
      (
       (metah (assoc 'meta description))
       (meta (cdr metah))
       (name (cdr (assoc 'name meta)))
       (ifaces (map lookup-interface (cdr (assoc 'interfaces meta))))
       (ifacedoc (map (lambda (iface)
			(ref-ify (symbol->string (car iface)))
			)
		      (reverse ifaces)))
       (engravers (filter
		   (lambda (x) (engraver-makes-grob? name x)) all-engravers-list))
       (namestr (symbol->string name))
       (engraver-names (map symbol->string (map ly:translator-name engravers)))
       )

    (make <texi-node>
      #:name namestr
      #:text
      (string-append
       namestr " grobs are created by: "
       (human-listify (map ref-ify
			   (map engraver-name engraver-names)))
       "\n\nStandard settings: \n\n"
       (grob-alist->texi description)
       "\n\nThis object supports the following interfaces: \n"
       (human-listify ifacedoc)
       ))
    ))

(define (all-grobs-doc)
  (make <texi-node>
    #:name "All layout objects"
    #:desc "Description and defaults for all Grobs"
    #:children
    (map (lambda (x) (grob-doc (cdr x)))  all-grob-descriptions)))

(define interface-description-alist
  (hash-fold
   (lambda (key val prior)
     (cons (cons key val)  prior)
     )
   '() (ly:all-grob-interfaces)))

(set! interface-description-alist (sort interface-description-alist alist<?))


;;;;;;;;;; check for dangling backend properties.
(define (mark-interface-properties entry)
  (map (lambda (x) (set-object-property! x  'iface-marked #t)) (caddr (cdr entry)))
  )

(map mark-interface-properties interface-description-alist)

(define (check-dangling-properties prop)
  (if (not (object-property prop 'iface-marked))
      (error  "\ndefine-grob-properties.scm: Can't find interface for property:" prop)))

(map check-dangling-properties all-backend-properties)

;;;;;;;;;;;;;;;;

(define (lookup-interface name)
  (let*  (
	  (entry  (hashq-ref (ly:all-grob-interfaces) name '() ))
	  )

    (if (equal? entry #f)
	(error "Unknown interface" name))
    
    entry
))

(define (all-interfaces-doc)
  (make <texi-node>
    #:name "Graphical Object Interfaces"
    #:desc "Building blocks of graphical objects"
    #:children
    (map interface-doc interface-description-alist)
    ))

(define (backend-properties-doc-string lst)
  (let*
      (
       (ps (sort (map symbol->string lst) string<?))
       (descs (map (lambda (prop)
		     (property->texi 'backend (string->symbol prop)  '()))
		   ps))
       (texi (description-list->texi descs))
       )
    texi))

  
;(dump-node (grob-doc (cdadr all-grob-descriptions))  (current-output-port) 0 )
(define (backend-doc-node)
  (make <texi-node>
    #:name "Backend"
    #:desc "Reference for the layout engine"
    #:children
    (list
     (all-grobs-doc)
     (all-interfaces-doc)
    (make <texi-node>
      #:name "User backend properties"
      #:desc "All tunable properties in a big list"
      #:text (backend-properties-doc-string all-internal-grob-properties))
    (make <texi-node>
      #:name "Internal backend properties"
      #:desc "All internal layout properties in a big list"
      #:text (backend-properties-doc-string all-user-grob-properties))
  )))
