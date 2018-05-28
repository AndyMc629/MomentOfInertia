from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.electronic_structure.plotter import DosPlotter

    
materials_list=["Cs3Pb2I7","Cs2RbPb2I7"]
for material in materials_list:
	dir = material + "/GeomOpt_" + material + "_I4mmm/BANDS/vasprun.xml"
	v = Vasprun(dir)
	bs = v.get_band_structure(line_mode=True)
	BSPlotter(bs).get_plot()
	BSPlotter(bs).save_plot(filename= material+"_band.eps")
      
	dosrun = Vasprun(material + "/GeomOpt_"+ material + "_I4mmm/DOS/vasprun.xml")
	dos = dosrun.complete_dos
	dp = DosPlotter(True, False,.2)
	delement = dos.get_element_dos()
	dp.add_dos_dict(delement)
	#dp.get_plot([-20,20]).show()
	dp.get_plot()
	dp.save_plot(filename=material+"_dos.eps")

dir="CsPbI3/GeomOpt_CsPbI3_pm3m/BANDS/vasprun.xml"
v = Vasprun(dir)
bs = v.get_band_structure(line_mode=True)
BSPlotter(bs).get_plot()
BSPlotter(bs).save_plot(filename= "CsPbI3_band.eps")

dosrun = Vasprun("CsPbI3/GeomOpt_CsPbI3_pm3m/DOS/vasprun.xml") 
dos = dosrun.complete_dos
dp = DosPlotter(True, False,.2)
delement = dos.get_element_dos()
dp.add_dos_dict(delement)
#dp.get_plot([-20,20]).show()
dp.get_plot()
dp.save_plot(filename="CsPbI3__dos.eps")
