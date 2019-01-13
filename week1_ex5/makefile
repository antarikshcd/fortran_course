#macro defining the target ie executable
target = runme.exe
# macro defininf the list of object files
# objects are recompiled in the order in which they appear in the list
# NOTE: TABS not required while defining macros
objs = m_data.o\
       m_alloc.o\
       main.o\
       s_initialize.o\
       s_fileout.o\
       s_updatefield.o\
       s_elem_updatefield.o\
       s_diagnostic.o\
       s_read_input.o

# linking: the target depends on the objects
# if objects are newer than the arget they are re-linked	
$(target):$(objs)
	f95 -ffree-form $(objs) -o $(target)

#dependencies:
main.o: main.f90 m_data.o m_alloc.o
	f95 -ffree-form -c main.f90
s_initialize.o: s_initialize.f90 m_data.o m_alloc.o
	f95 -ffree-form -c s_initialize.f90
s_fileout.o: s_fileout.f90
	f95 -ffree-form -c s_fileout.f90
s_updatefield.o: s_updatefield.f90
	f95 -ffree-form -c s_updatefield.f90
s_elem_updatefield.o: s_elem_updatefield.f90
	f95 -ffree-form -c s_elem_updatefield.f90	
s_diagnostic.o: s_diagnostic.f90
	f95 -ffree-form -c s_diagnostic.f90
s_read_input.o: s_read_input.f90
	f95 -ffree-form -c s_read_input.f90	
m_alloc.o: m_alloc.f90 m_data.o
	f95 -ffree-form -c m_alloc.f90
m_data.o: m_data.f90
	f95 -ffree-form -c m_data.f90

new: clean $(target)
clean:
	rm -fr $(objs)
	rm -fr *.dat		

