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
rz(-2.2890685) q[0];
sx q[0];
rz(-0.26242119) q[0];
rz(2.8334795) q[1];
sx q[1];
rz(-2.3051655) q[1];
sx q[1];
rz(3.1321373) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1605837) q[0];
sx q[0];
rz(-2.0045337) q[0];
sx q[0];
rz(-1.0354614) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1310102) q[2];
sx q[2];
rz(-0.081801266) q[2];
sx q[2];
rz(0.33120016) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7345924) q[1];
sx q[1];
rz(-1.0838036) q[1];
sx q[1];
rz(-0.89501801) q[1];
rz(-0.81266021) q[3];
sx q[3];
rz(-1.9403807) q[3];
sx q[3];
rz(0.11358914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1001763) q[2];
sx q[2];
rz(-0.88465038) q[2];
sx q[2];
rz(0.46119383) q[2];
rz(2.6395116) q[3];
sx q[3];
rz(-0.82379782) q[3];
sx q[3];
rz(1.4095149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0111888) q[0];
sx q[0];
rz(-1.4606425) q[0];
sx q[0];
rz(0.23656626) q[0];
rz(0.81831167) q[1];
sx q[1];
rz(-1.9383483) q[1];
sx q[1];
rz(2.1990105) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.062898) q[0];
sx q[0];
rz(-1.601614) q[0];
sx q[0];
rz(0.0027058733) q[0];
x q[1];
rz(1.0636395) q[2];
sx q[2];
rz(-0.13826577) q[2];
sx q[2];
rz(-2.0035715) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.67737867) q[1];
sx q[1];
rz(-0.83589743) q[1];
sx q[1];
rz(-1.5688495) q[1];
rz(-0.97901042) q[3];
sx q[3];
rz(-0.56603449) q[3];
sx q[3];
rz(-2.4477959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73611034) q[2];
sx q[2];
rz(-0.63850275) q[2];
sx q[2];
rz(0.68518266) q[2];
rz(-2.3403366) q[3];
sx q[3];
rz(-1.6359436) q[3];
sx q[3];
rz(1.8906458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7291229) q[0];
sx q[0];
rz(-0.47054371) q[0];
sx q[0];
rz(-2.1614918) q[0];
rz(0.12067548) q[1];
sx q[1];
rz(-1.1680892) q[1];
sx q[1];
rz(-1.6623704) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2187933) q[0];
sx q[0];
rz(-2.3258219) q[0];
sx q[0];
rz(1.3479309) q[0];
rz(-pi) q[1];
rz(-2.4083869) q[2];
sx q[2];
rz(-0.63989675) q[2];
sx q[2];
rz(-2.899183) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4795717) q[1];
sx q[1];
rz(-1.5618854) q[1];
sx q[1];
rz(-1.5655976) q[1];
x q[2];
rz(0.48101403) q[3];
sx q[3];
rz(-2.4949007) q[3];
sx q[3];
rz(-2.8493136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4270619) q[2];
sx q[2];
rz(-1.6468628) q[2];
sx q[2];
rz(-2.96116) q[2];
rz(2.7212972) q[3];
sx q[3];
rz(-2.1642978) q[3];
sx q[3];
rz(2.596627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20035289) q[0];
sx q[0];
rz(-1.350133) q[0];
sx q[0];
rz(-3.0923162) q[0];
rz(2.409528) q[1];
sx q[1];
rz(-2.2923636) q[1];
sx q[1];
rz(-1.3759618) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7017562) q[0];
sx q[0];
rz(-1.4539833) q[0];
sx q[0];
rz(-0.18420903) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76566546) q[2];
sx q[2];
rz(-1.7244974) q[2];
sx q[2];
rz(-1.1534899) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1092618) q[1];
sx q[1];
rz(-1.9606661) q[1];
sx q[1];
rz(1.6435739) q[1];
x q[2];
rz(-0.64582156) q[3];
sx q[3];
rz(-2.1649515) q[3];
sx q[3];
rz(2.770383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.10776082) q[2];
sx q[2];
rz(-1.3904479) q[2];
sx q[2];
rz(-2.4793009) q[2];
rz(-1.8108588) q[3];
sx q[3];
rz(-0.8225421) q[3];
sx q[3];
rz(1.2983769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0419643) q[0];
sx q[0];
rz(-2.8350416) q[0];
sx q[0];
rz(-2.1339259) q[0];
rz(3.0362466) q[1];
sx q[1];
rz(-2.4844929) q[1];
sx q[1];
rz(2.9109921) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14734257) q[0];
sx q[0];
rz(-1.2201637) q[0];
sx q[0];
rz(-2.7870534) q[0];
x q[1];
rz(1.2097539) q[2];
sx q[2];
rz(-1.3461543) q[2];
sx q[2];
rz(1.063907) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.19781348) q[1];
sx q[1];
rz(-2.9557485) q[1];
sx q[1];
rz(0.4798823) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6449987) q[3];
sx q[3];
rz(-1.4497744) q[3];
sx q[3];
rz(-1.6337486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0838123) q[2];
sx q[2];
rz(-0.96692204) q[2];
sx q[2];
rz(-0.18873611) q[2];
rz(1.905929) q[3];
sx q[3];
rz(-2.6344447) q[3];
sx q[3];
rz(-1.2602826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96596232) q[0];
sx q[0];
rz(-1.2164793) q[0];
sx q[0];
rz(2.8438582) q[0];
rz(3.0689902) q[1];
sx q[1];
rz(-0.68521348) q[1];
sx q[1];
rz(2.4593478) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84083623) q[0];
sx q[0];
rz(-1.8634691) q[0];
sx q[0];
rz(0.93772823) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3700245) q[2];
sx q[2];
rz(-2.4983495) q[2];
sx q[2];
rz(2.7386576) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1156716) q[1];
sx q[1];
rz(-2.762315) q[1];
sx q[1];
rz(1.161452) q[1];
rz(-pi) q[2];
rz(-2.0892049) q[3];
sx q[3];
rz(-2.4757407) q[3];
sx q[3];
rz(-0.4308373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5630774) q[2];
sx q[2];
rz(-0.35334057) q[2];
sx q[2];
rz(-3.0541218) q[2];
rz(0.89720094) q[3];
sx q[3];
rz(-1.2465435) q[3];
sx q[3];
rz(-2.7520666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099365756) q[0];
sx q[0];
rz(-1.1023738) q[0];
sx q[0];
rz(1.9788096) q[0];
rz(0.21576628) q[1];
sx q[1];
rz(-0.81753221) q[1];
sx q[1];
rz(-2.20631) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29998818) q[0];
sx q[0];
rz(-0.81735742) q[0];
sx q[0];
rz(1.2410844) q[0];
x q[1];
rz(-1.0583508) q[2];
sx q[2];
rz(-2.847306) q[2];
sx q[2];
rz(2.2541915) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.062309655) q[1];
sx q[1];
rz(-2.3543962) q[1];
sx q[1];
rz(1.4592378) q[1];
rz(0.8029253) q[3];
sx q[3];
rz(-0.438779) q[3];
sx q[3];
rz(-0.9641274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0995348) q[2];
sx q[2];
rz(-2.4877986) q[2];
sx q[2];
rz(-2.330244) q[2];
rz(2.8387496) q[3];
sx q[3];
rz(-1.1648014) q[3];
sx q[3];
rz(2.5109049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52043668) q[0];
sx q[0];
rz(-2.8148837) q[0];
sx q[0];
rz(-0.69123554) q[0];
rz(0.65903819) q[1];
sx q[1];
rz(-1.5182511) q[1];
sx q[1];
rz(2.5504327) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49159971) q[0];
sx q[0];
rz(-0.2831471) q[0];
sx q[0];
rz(1.6121907) q[0];
rz(0.28556683) q[2];
sx q[2];
rz(-1.0896249) q[2];
sx q[2];
rz(-2.8701631) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8670676) q[1];
sx q[1];
rz(-2.3127416) q[1];
sx q[1];
rz(-1.1222003) q[1];
x q[2];
rz(1.2783613) q[3];
sx q[3];
rz(-1.5163444) q[3];
sx q[3];
rz(-1.4002348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.14472321) q[2];
sx q[2];
rz(-2.150712) q[2];
sx q[2];
rz(1.4608176) q[2];
rz(2.3878494) q[3];
sx q[3];
rz(-1.6921348) q[3];
sx q[3];
rz(-0.47372216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4562456) q[0];
sx q[0];
rz(-2.0401177) q[0];
sx q[0];
rz(-2.9154678) q[0];
rz(1.6926758) q[1];
sx q[1];
rz(-2.1782404) q[1];
sx q[1];
rz(-1.0015782) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6152214) q[0];
sx q[0];
rz(-0.60524056) q[0];
sx q[0];
rz(0.25608082) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4422574) q[2];
sx q[2];
rz(-1.275389) q[2];
sx q[2];
rz(0.70198529) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.68943095) q[1];
sx q[1];
rz(-0.2161444) q[1];
sx q[1];
rz(-0.5139022) q[1];
x q[2];
rz(-0.36205451) q[3];
sx q[3];
rz(-2.4241862) q[3];
sx q[3];
rz(-0.91594726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.61233026) q[2];
sx q[2];
rz(-2.3713106) q[2];
sx q[2];
rz(2.9840577) q[2];
rz(0.15466386) q[3];
sx q[3];
rz(-1.9779343) q[3];
sx q[3];
rz(0.4304339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4100819) q[0];
sx q[0];
rz(-1.9022576) q[0];
sx q[0];
rz(2.4391158) q[0];
rz(-0.53600535) q[1];
sx q[1];
rz(-2.3119226) q[1];
sx q[1];
rz(-1.3943256) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0939801) q[0];
sx q[0];
rz(-1.5965921) q[0];
sx q[0];
rz(1.4071147) q[0];
rz(-pi) q[1];
x q[1];
rz(1.461566) q[2];
sx q[2];
rz(-1.6642206) q[2];
sx q[2];
rz(2.1364853) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.9623757) q[1];
sx q[1];
rz(-2.0834647) q[1];
sx q[1];
rz(-2.8151399) q[1];
rz(-pi) q[2];
rz(1.5709747) q[3];
sx q[3];
rz(-1.1552253) q[3];
sx q[3];
rz(3.0316169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.006762) q[2];
sx q[2];
rz(-2.4595021) q[2];
sx q[2];
rz(1.6949867) q[2];
rz(-0.26099482) q[3];
sx q[3];
rz(-2.1311396) q[3];
sx q[3];
rz(-1.9058156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9217054) q[0];
sx q[0];
rz(-0.6548665) q[0];
sx q[0];
rz(-1.0153216) q[0];
rz(-2.0657397) q[1];
sx q[1];
rz(-1.8668108) q[1];
sx q[1];
rz(1.6783953) q[1];
rz(2.867472) q[2];
sx q[2];
rz(-2.7514396) q[2];
sx q[2];
rz(-0.20382602) q[2];
rz(0.76293972) q[3];
sx q[3];
rz(-1.2423824) q[3];
sx q[3];
rz(0.35425274) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
