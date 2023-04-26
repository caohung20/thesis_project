import sparse_matrix_lfnewton as lfnewton
tech=lfnewton.tech
nbus=lfnewton.nbus
deltad=lfnewton.deltad
Vm=lfnewton.Vm
Pd=lfnewton.Pd
Qd=lfnewton.Qd
Pg=lfnewton.Pg
Qg=lfnewton.Qg
Qsh=lfnewton.Qsh
Pdt=lfnewton.Pdt
Qdt=lfnewton.Qdt
Pgt=lfnewton.Pgt
Qgt=lfnewton.Qgt
Qsht=lfnewton.Qsht
maxerror=lfnewton.maxerror
iter=lfnewton.iter
print(tech)
print('                      Maximum Power Mismatch = %g \n' % maxerror)
print('                             No. of Iterations = %g \n\n' % iter)

head = [    '    Bus  Voltage  Angle    ------Load------    ---Generation---   Injected',    '    No.  Mag.     Degree     MW       Mvar       MW       Mvar       Mvar ',    '                                                                          ']
print('\n'.join(head))

for n in range(nbus):
    print(' %5g' % n, end=' ')
    print(' %7.3f' % Vm[n].real, end=' ')
    print(' %8.3f' % deltad[n].real, end=' ')
    print(' %9.3f' % Pd[n], end=' ')
    print(' %9.3f' % Qd[n], end=' ')
    print(' %9.3f' % Pg[n].real, end=' ')
    print(' %9.3f ' % Qg[n].real, end=' ')
    print(' %8.3f' % Qsh[n])

print('      ')
print('    Total              ', end=' ')
print(' %9.3f' % Pdt, end=' ')
print(' %9.3f' % Qdt, end=' ')
print(' %9.3f' % Pgt.real, end=' ')
print(' %9.3f' % Qgt.real, end=' ')
print(' %9.3f\n\n' % Qsht)
