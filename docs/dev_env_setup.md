Setup development environment
=============================

Date: Jul 23, 2015

  1. Install [anaconda](http://continuum.io/downloads)

  2. Create __imodels__ environment:
  
  		$ conda create -n imodels python=3.4
  		
  3. Activate __imodels__ environment:

        $ source activate imodels
        
  4. Install some conda packages:
  
        $ conda install pkgconfig matplotlib networkx ipython ipython-notebook
        $ conda install -c omnia mdtraj
        
  5. Patch _mdtraj_ with [this patch](https://github.com/mdtraj/mdtraj/pull/896/files).
  
  		$ cd /path/to/anaconda/envs/imodels/lib/python3.4/site-packages/mdtraj/html
  		$ subl trajectory_widget.py
  		
  6. Clone [Intermol](https://github.com/shirtsgroup/InterMol), [foyer](https://github.com/iModels/foyer) and [mbuild](https://github.com/iModels/mbuild) projects and setup them:
  
  		$ cd /path/to/your/projects/
  		$ git clone https://github.com/shirtsgroup/InterMol.git
  		$ git clone https://github.com/iModels/foyer.git
  		$ git clone https://github.com/iModels/mbuild.git

	##### Intermol
  		$ cd Intermol
  		$ git checkout forcetype
  		$ python setup.py develop
  		
  	##### foyer
  		$ cd foyer
  		$ python setup.py develop
  		
  	##### mbuild
  		$ cd mbuild
  		$ python setup.py develop
  	
  	