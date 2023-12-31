OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7744301) q[0];
sx q[0];
rz(-0.91355938) q[0];
sx q[0];
rz(1.4120742) q[0];
rz(-2.9867759) q[1];
sx q[1];
rz(-0.59564367) q[1];
sx q[1];
rz(-1.6593978) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5858551) q[0];
sx q[0];
rz(-1.9940388) q[0];
sx q[0];
rz(-1.7413571) q[0];
x q[1];
rz(2.8843845) q[2];
sx q[2];
rz(-1.6991985) q[2];
sx q[2];
rz(-0.53127015) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2130148) q[1];
sx q[1];
rz(-1.9933812) q[1];
sx q[1];
rz(-2.1133452) q[1];
x q[2];
rz(-1.0079908) q[3];
sx q[3];
rz(-1.4853012) q[3];
sx q[3];
rz(0.29719719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1564864) q[2];
sx q[2];
rz(-0.50922314) q[2];
sx q[2];
rz(0.86581725) q[2];
rz(-2.1872897) q[3];
sx q[3];
rz(-1.6031957) q[3];
sx q[3];
rz(1.2877864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99825478) q[0];
sx q[0];
rz(-1.4366432) q[0];
sx q[0];
rz(-3.1153733) q[0];
rz(-1.5401309) q[1];
sx q[1];
rz(-1.5427579) q[1];
sx q[1];
rz(2.1781133) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.119495) q[0];
sx q[0];
rz(-2.1848328) q[0];
sx q[0];
rz(-0.0028145785) q[0];
rz(-2.8112667) q[2];
sx q[2];
rz(-2.1147554) q[2];
sx q[2];
rz(-0.75737539) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.023167921) q[1];
sx q[1];
rz(-1.8903362) q[1];
sx q[1];
rz(-2.097514) q[1];
rz(-0.96879949) q[3];
sx q[3];
rz(-0.41799823) q[3];
sx q[3];
rz(0.19817643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6271237) q[2];
sx q[2];
rz(-1.1273948) q[2];
sx q[2];
rz(0.13452402) q[2];
rz(2.3965805) q[3];
sx q[3];
rz(-2.9146505) q[3];
sx q[3];
rz(-2.1988595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2117675) q[0];
sx q[0];
rz(-2.7524502) q[0];
sx q[0];
rz(2.3441558) q[0];
rz(-2.0939317) q[1];
sx q[1];
rz(-0.14973775) q[1];
sx q[1];
rz(-2.581596) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5138429) q[0];
sx q[0];
rz(-1.3686413) q[0];
sx q[0];
rz(1.1802799) q[0];
rz(-0.067990818) q[2];
sx q[2];
rz(-2.8376841) q[2];
sx q[2];
rz(-1.4246724) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8078976) q[1];
sx q[1];
rz(-1.9296608) q[1];
sx q[1];
rz(1.3586504) q[1];
x q[2];
rz(1.0873763) q[3];
sx q[3];
rz(-1.1421575) q[3];
sx q[3];
rz(-2.7408858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3893163) q[2];
sx q[2];
rz(-1.1976778) q[2];
sx q[2];
rz(2.9690572) q[2];
rz(-2.1595188) q[3];
sx q[3];
rz(-1.3970102) q[3];
sx q[3];
rz(1.0579695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.003222) q[0];
sx q[0];
rz(-1.0976185) q[0];
sx q[0];
rz(-0.28451434) q[0];
rz(2.8248887) q[1];
sx q[1];
rz(-0.43276325) q[1];
sx q[1];
rz(1.8428615) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1176778) q[0];
sx q[0];
rz(-1.0751343) q[0];
sx q[0];
rz(0.25650521) q[0];
x q[1];
rz(-1.5721333) q[2];
sx q[2];
rz(-0.7664116) q[2];
sx q[2];
rz(3.1198451) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.58060927) q[1];
sx q[1];
rz(-2.1316075) q[1];
sx q[1];
rz(2.9210655) q[1];
x q[2];
rz(3.0292547) q[3];
sx q[3];
rz(-1.897052) q[3];
sx q[3];
rz(-2.5374967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6056885) q[2];
sx q[2];
rz(-0.19583344) q[2];
sx q[2];
rz(-2.7569125) q[2];
rz(0.7540594) q[3];
sx q[3];
rz(-1.056517) q[3];
sx q[3];
rz(1.6872905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5383179) q[0];
sx q[0];
rz(-0.98709995) q[0];
sx q[0];
rz(-1.3866562) q[0];
rz(0.23100135) q[1];
sx q[1];
rz(-1.8004386) q[1];
sx q[1];
rz(2.8447661) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0083858) q[0];
sx q[0];
rz(-1.0002245) q[0];
sx q[0];
rz(1.5242566) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9752521) q[2];
sx q[2];
rz(-2.3918249) q[2];
sx q[2];
rz(-1.9094085) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0592812) q[1];
sx q[1];
rz(-1.1059522) q[1];
sx q[1];
rz(2.6962198) q[1];
rz(-pi) q[2];
rz(0.71367587) q[3];
sx q[3];
rz(-1.6380966) q[3];
sx q[3];
rz(2.0590559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37830535) q[2];
sx q[2];
rz(-1.3102691) q[2];
sx q[2];
rz(-2.7491167) q[2];
rz(1.1522419) q[3];
sx q[3];
rz(-2.4270054) q[3];
sx q[3];
rz(2.8241482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0157938) q[0];
sx q[0];
rz(-1.5725461) q[0];
sx q[0];
rz(-2.3902067) q[0];
rz(-1.8136576) q[1];
sx q[1];
rz(-1.8782047) q[1];
sx q[1];
rz(2.5352535) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060121814) q[0];
sx q[0];
rz(-1.524964) q[0];
sx q[0];
rz(1.5874552) q[0];
rz(-pi) q[1];
rz(-2.0727856) q[2];
sx q[2];
rz(-0.74337465) q[2];
sx q[2];
rz(2.5943287) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1746754) q[1];
sx q[1];
rz(-1.3904018) q[1];
sx q[1];
rz(-2.128152) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0665928) q[3];
sx q[3];
rz(-1.3538176) q[3];
sx q[3];
rz(-1.0523588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5027344) q[2];
sx q[2];
rz(-2.0998462) q[2];
sx q[2];
rz(2.0022557) q[2];
rz(1.4849439) q[3];
sx q[3];
rz(-1.9610201) q[3];
sx q[3];
rz(-3.0373354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5320324) q[0];
sx q[0];
rz(-2.39344) q[0];
sx q[0];
rz(0.50810057) q[0];
rz(1.5787026) q[1];
sx q[1];
rz(-2.0527614) q[1];
sx q[1];
rz(-0.79024822) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3797487) q[0];
sx q[0];
rz(-1.8659235) q[0];
sx q[0];
rz(-0.6167114) q[0];
rz(1.5590645) q[2];
sx q[2];
rz(-1.7844229) q[2];
sx q[2];
rz(-0.27660433) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33752791) q[1];
sx q[1];
rz(-1.3009326) q[1];
sx q[1];
rz(-2.7359664) q[1];
rz(-pi) q[2];
rz(0.64738691) q[3];
sx q[3];
rz(-1.2643407) q[3];
sx q[3];
rz(-2.2341773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.053085176) q[2];
sx q[2];
rz(-2.6997456) q[2];
sx q[2];
rz(-1.7283758) q[2];
rz(-1.6648071) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(-3.0800381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4247894) q[0];
sx q[0];
rz(-3.1112818) q[0];
sx q[0];
rz(1.0472263) q[0];
rz(-0.60910243) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(-1.3887127) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0213288) q[0];
sx q[0];
rz(-1.6006032) q[0];
sx q[0];
rz(2.1304312) q[0];
rz(-pi) q[1];
rz(-2.6128204) q[2];
sx q[2];
rz(-1.1985589) q[2];
sx q[2];
rz(-2.1793274) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.82909225) q[1];
sx q[1];
rz(-1.478985) q[1];
sx q[1];
rz(-1.5822259) q[1];
rz(-pi) q[2];
rz(-2.1771168) q[3];
sx q[3];
rz(-1.3250321) q[3];
sx q[3];
rz(2.5442459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1887112) q[2];
sx q[2];
rz(-0.41023508) q[2];
sx q[2];
rz(-0.88225538) q[2];
rz(-1.7404209) q[3];
sx q[3];
rz(-1.975235) q[3];
sx q[3];
rz(1.2004948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8700478) q[0];
sx q[0];
rz(-2.7247868) q[0];
sx q[0];
rz(-1.7154988) q[0];
rz(3.0601314) q[1];
sx q[1];
rz(-1.9790244) q[1];
sx q[1];
rz(-2.5833599) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4687913) q[0];
sx q[0];
rz(-1.1835915) q[0];
sx q[0];
rz(-1.1084491) q[0];
rz(-pi) q[1];
rz(1.105537) q[2];
sx q[2];
rz(-0.12013398) q[2];
sx q[2];
rz(2.5533954) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7575175) q[1];
sx q[1];
rz(-1.6346524) q[1];
sx q[1];
rz(-0.82660316) q[1];
x q[2];
rz(2.8006323) q[3];
sx q[3];
rz(-2.6302359) q[3];
sx q[3];
rz(-2.2380059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2925064) q[2];
sx q[2];
rz(-1.2693274) q[2];
sx q[2];
rz(-0.212184) q[2];
rz(0.21197453) q[3];
sx q[3];
rz(-0.68325716) q[3];
sx q[3];
rz(1.2020483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6367209) q[0];
sx q[0];
rz(-2.3175406) q[0];
sx q[0];
rz(1.5378392) q[0];
rz(0.82540712) q[1];
sx q[1];
rz(-2.4688265) q[1];
sx q[1];
rz(-2.6182981) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035318035) q[0];
sx q[0];
rz(-0.9629074) q[0];
sx q[0];
rz(1.4950698) q[0];
rz(0.9988437) q[2];
sx q[2];
rz(-1.5339282) q[2];
sx q[2];
rz(1.8429304) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.938898) q[1];
sx q[1];
rz(-1.2564066) q[1];
sx q[1];
rz(-0.37757341) q[1];
x q[2];
rz(1.80199) q[3];
sx q[3];
rz(-1.8095067) q[3];
sx q[3];
rz(1.064144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3045197) q[2];
sx q[2];
rz(-0.7154811) q[2];
sx q[2];
rz(0.26930299) q[2];
rz(2.6473911) q[3];
sx q[3];
rz(-2.2952357) q[3];
sx q[3];
rz(2.0555029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3257278) q[0];
sx q[0];
rz(-1.5300735) q[0];
sx q[0];
rz(-1.6515401) q[0];
rz(1.6745463) q[1];
sx q[1];
rz(-0.29232262) q[1];
sx q[1];
rz(1.2437337) q[1];
rz(-1.5291443) q[2];
sx q[2];
rz(-0.26298444) q[2];
sx q[2];
rz(2.0900805) q[2];
rz(-2.9071517) q[3];
sx q[3];
rz(-0.81080484) q[3];
sx q[3];
rz(-0.57639359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
