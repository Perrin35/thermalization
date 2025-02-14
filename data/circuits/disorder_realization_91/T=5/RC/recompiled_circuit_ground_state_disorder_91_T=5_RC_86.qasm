OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.087091669) q[0];
sx q[0];
rz(0.85252419) q[0];
sx q[0];
rz(9.6871992) q[0];
rz(2.8334795) q[1];
sx q[1];
rz(-2.3051655) q[1];
sx q[1];
rz(-0.0094553789) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1605837) q[0];
sx q[0];
rz(-1.137059) q[0];
sx q[0];
rz(2.1061312) q[0];
rz(-pi) q[1];
rz(-1.6401372) q[2];
sx q[2];
rz(-1.5273646) q[2];
sx q[2];
rz(-0.68088898) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5227185) q[1];
sx q[1];
rz(-2.1565003) q[1];
sx q[1];
rz(2.5453091) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3289324) q[3];
sx q[3];
rz(-1.201212) q[3];
sx q[3];
rz(-0.11358914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0414163) q[2];
sx q[2];
rz(-0.88465038) q[2];
sx q[2];
rz(0.46119383) q[2];
rz(0.50208107) q[3];
sx q[3];
rz(-2.3177948) q[3];
sx q[3];
rz(1.4095149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0111888) q[0];
sx q[0];
rz(-1.6809502) q[0];
sx q[0];
rz(0.23656626) q[0];
rz(-2.323281) q[1];
sx q[1];
rz(-1.9383483) q[1];
sx q[1];
rz(2.1990105) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1326951) q[0];
sx q[0];
rz(-3.1106565) q[0];
sx q[0];
rz(-1.6583468) q[0];
x q[1];
rz(-0.06748345) q[2];
sx q[2];
rz(-1.6915671) q[2];
sx q[2];
rz(-1.4923332) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.68028223) q[1];
sx q[1];
rz(-0.73490099) q[1];
sx q[1];
rz(-0.0021541455) q[1];
rz(-pi) q[2];
rz(-2.8009612) q[3];
sx q[3];
rz(-1.1095205) q[3];
sx q[3];
rz(-0.021286437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4054823) q[2];
sx q[2];
rz(-0.63850275) q[2];
sx q[2];
rz(-0.68518266) q[2];
rz(-0.80125609) q[3];
sx q[3];
rz(-1.5056491) q[3];
sx q[3];
rz(1.8906458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41246978) q[0];
sx q[0];
rz(-0.47054371) q[0];
sx q[0];
rz(2.1614918) q[0];
rz(-3.0209172) q[1];
sx q[1];
rz(-1.1680892) q[1];
sx q[1];
rz(-1.6623704) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89941601) q[0];
sx q[0];
rz(-0.78095312) q[0];
sx q[0];
rz(-2.910898) q[0];
x q[1];
rz(0.50522106) q[2];
sx q[2];
rz(-1.9818993) q[2];
sx q[2];
rz(1.9540602) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66202098) q[1];
sx q[1];
rz(-1.5797073) q[1];
sx q[1];
rz(1.5655976) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2347264) q[3];
sx q[3];
rz(-2.1343291) q[3];
sx q[3];
rz(0.8714217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4270619) q[2];
sx q[2];
rz(-1.4947299) q[2];
sx q[2];
rz(-2.96116) q[2];
rz(2.7212972) q[3];
sx q[3];
rz(-2.1642978) q[3];
sx q[3];
rz(-0.54496566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9412398) q[0];
sx q[0];
rz(-1.350133) q[0];
sx q[0];
rz(-0.049276503) q[0];
rz(-0.73206466) q[1];
sx q[1];
rz(-0.8492291) q[1];
sx q[1];
rz(-1.7656309) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7017562) q[0];
sx q[0];
rz(-1.4539833) q[0];
sx q[0];
rz(2.9573836) q[0];
rz(2.3759272) q[2];
sx q[2];
rz(-1.7244974) q[2];
sx q[2];
rz(-1.1534899) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84280864) q[1];
sx q[1];
rz(-2.7453303) q[1];
sx q[1];
rz(2.9664459) q[1];
rz(-2.4957711) q[3];
sx q[3];
rz(-0.97664112) q[3];
sx q[3];
rz(2.770383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0338318) q[2];
sx q[2];
rz(-1.3904479) q[2];
sx q[2];
rz(2.4793009) q[2];
rz(-1.8108588) q[3];
sx q[3];
rz(-2.3190506) q[3];
sx q[3];
rz(1.8432157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0419643) q[0];
sx q[0];
rz(-0.30655107) q[0];
sx q[0];
rz(1.0076667) q[0];
rz(-0.10534605) q[1];
sx q[1];
rz(-2.4844929) q[1];
sx q[1];
rz(-0.23060051) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67543519) q[0];
sx q[0];
rz(-2.648232) q[0];
sx q[0];
rz(-0.81146474) q[0];
x q[1];
rz(1.2097539) q[2];
sx q[2];
rz(-1.3461543) q[2];
sx q[2];
rz(1.063907) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.28921738) q[1];
sx q[1];
rz(-1.7354449) q[1];
sx q[1];
rz(1.4842121) q[1];
x q[2];
rz(2.9419961) q[3];
sx q[3];
rz(-2.4869339) q[3];
sx q[3];
rz(-2.9195291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0577804) q[2];
sx q[2];
rz(-0.96692204) q[2];
sx q[2];
rz(-2.9528565) q[2];
rz(-1.905929) q[3];
sx q[3];
rz(-0.50714791) q[3];
sx q[3];
rz(-1.2602826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96596232) q[0];
sx q[0];
rz(-1.9251134) q[0];
sx q[0];
rz(-0.29773444) q[0];
rz(-3.0689902) q[1];
sx q[1];
rz(-0.68521348) q[1];
sx q[1];
rz(0.6822449) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1045584) q[0];
sx q[0];
rz(-0.68892043) q[0];
sx q[0];
rz(1.0997186) q[0];
rz(-0.14839006) q[2];
sx q[2];
rz(-0.94252693) q[2];
sx q[2];
rz(2.4895957) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9798237) q[1];
sx q[1];
rz(-1.7186972) q[1];
sx q[1];
rz(-1.921341) q[1];
rz(0.37118427) q[3];
sx q[3];
rz(-2.1371578) q[3];
sx q[3];
rz(2.0831828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5630774) q[2];
sx q[2];
rz(-0.35334057) q[2];
sx q[2];
rz(-3.0541218) q[2];
rz(-2.2443917) q[3];
sx q[3];
rz(-1.8950491) q[3];
sx q[3];
rz(2.7520666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0422269) q[0];
sx q[0];
rz(-1.1023738) q[0];
sx q[0];
rz(1.1627831) q[0];
rz(2.9258264) q[1];
sx q[1];
rz(-0.81753221) q[1];
sx q[1];
rz(2.20631) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29998818) q[0];
sx q[0];
rz(-0.81735742) q[0];
sx q[0];
rz(1.2410844) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0832418) q[2];
sx q[2];
rz(-2.847306) q[2];
sx q[2];
rz(0.8874011) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.079283) q[1];
sx q[1];
rz(-0.78719646) q[1];
sx q[1];
rz(-1.4592378) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8264826) q[3];
sx q[3];
rz(-1.8813881) q[3];
sx q[3];
rz(0.14652182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0995348) q[2];
sx q[2];
rz(-2.4877986) q[2];
sx q[2];
rz(-0.81134861) q[2];
rz(0.30284303) q[3];
sx q[3];
rz(-1.1648014) q[3];
sx q[3];
rz(-2.5109049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52043668) q[0];
sx q[0];
rz(-2.8148837) q[0];
sx q[0];
rz(2.4503571) q[0];
rz(0.65903819) q[1];
sx q[1];
rz(-1.5182511) q[1];
sx q[1];
rz(-0.59115994) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6499929) q[0];
sx q[0];
rz(-0.2831471) q[0];
sx q[0];
rz(1.6121907) q[0];
rz(-pi) q[1];
rz(2.8560258) q[2];
sx q[2];
rz(-1.0896249) q[2];
sx q[2];
rz(-0.27142957) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2745251) q[1];
sx q[1];
rz(-0.82885107) q[1];
sx q[1];
rz(-1.1222003) q[1];
x q[2];
rz(0.05686111) q[3];
sx q[3];
rz(-1.8627852) q[3];
sx q[3];
rz(-0.1869456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.14472321) q[2];
sx q[2];
rz(-0.99088061) q[2];
sx q[2];
rz(-1.4608176) q[2];
rz(-2.3878494) q[3];
sx q[3];
rz(-1.6921348) q[3];
sx q[3];
rz(-2.6678705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4562456) q[0];
sx q[0];
rz(-2.0401177) q[0];
sx q[0];
rz(0.22612485) q[0];
rz(1.6926758) q[1];
sx q[1];
rz(-0.96335226) q[1];
sx q[1];
rz(1.0015782) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9739363) q[0];
sx q[0];
rz(-1.7154124) q[0];
sx q[0];
rz(-2.5517716) q[0];
rz(1.1923157) q[2];
sx q[2];
rz(-0.90736249) q[2];
sx q[2];
rz(-1.1089693) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4521617) q[1];
sx q[1];
rz(-2.9254483) q[1];
sx q[1];
rz(-2.6276905) q[1];
x q[2];
rz(-2.4572152) q[3];
sx q[3];
rz(-1.805814) q[3];
sx q[3];
rz(-0.37684611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5292624) q[2];
sx q[2];
rz(-2.3713106) q[2];
sx q[2];
rz(-2.9840577) q[2];
rz(-0.15466386) q[3];
sx q[3];
rz(-1.1636584) q[3];
sx q[3];
rz(0.4304339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73151076) q[0];
sx q[0];
rz(-1.2393351) q[0];
sx q[0];
rz(-2.4391158) q[0];
rz(2.6055873) q[1];
sx q[1];
rz(-0.82967007) q[1];
sx q[1];
rz(-1.747267) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6226688) q[0];
sx q[0];
rz(-1.4071696) q[0];
sx q[0];
rz(-0.026145025) q[0];
rz(-0.86081204) q[2];
sx q[2];
rz(-0.14361266) q[2];
sx q[2];
rz(-1.8710305) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.9623757) q[1];
sx q[1];
rz(-2.0834647) q[1];
sx q[1];
rz(-2.8151399) q[1];
rz(-3.1411885) q[3];
sx q[3];
rz(-2.7260216) q[3];
sx q[3];
rz(3.0311751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.006762) q[2];
sx q[2];
rz(-2.4595021) q[2];
sx q[2];
rz(-1.4466059) q[2];
rz(-0.26099482) q[3];
sx q[3];
rz(-2.1311396) q[3];
sx q[3];
rz(1.235777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9217054) q[0];
sx q[0];
rz(-0.6548665) q[0];
sx q[0];
rz(-1.0153216) q[0];
rz(1.075853) q[1];
sx q[1];
rz(-1.8668108) q[1];
sx q[1];
rz(1.6783953) q[1];
rz(0.27412065) q[2];
sx q[2];
rz(-0.39015301) q[2];
sx q[2];
rz(2.9377666) q[2];
rz(-2.0113284) q[3];
sx q[3];
rz(-0.85778271) q[3];
sx q[3];
rz(-1.5156619) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
