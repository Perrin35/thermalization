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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.981009) q[0];
sx q[0];
rz(-1.137059) q[0];
sx q[0];
rz(1.0354614) q[0];
rz(-1.6401372) q[2];
sx q[2];
rz(-1.614228) q[2];
sx q[2];
rz(0.68088898) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7345924) q[1];
sx q[1];
rz(-2.0577891) q[1];
sx q[1];
rz(0.89501801) q[1];
rz(2.0838787) q[3];
sx q[3];
rz(-2.3145967) q[3];
sx q[3];
rz(1.3199453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0414163) q[2];
sx q[2];
rz(-0.88465038) q[2];
sx q[2];
rz(-0.46119383) q[2];
rz(2.6395116) q[3];
sx q[3];
rz(-2.3177948) q[3];
sx q[3];
rz(-1.4095149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0111888) q[0];
sx q[0];
rz(-1.6809502) q[0];
sx q[0];
rz(2.9050264) q[0];
rz(0.81831167) q[1];
sx q[1];
rz(-1.2032443) q[1];
sx q[1];
rz(0.94258211) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0088975178) q[0];
sx q[0];
rz(-0.030936154) q[0];
sx q[0];
rz(1.6583468) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0741092) q[2];
sx q[2];
rz(-1.4500256) q[2];
sx q[2];
rz(1.4923332) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.68028223) q[1];
sx q[1];
rz(-0.73490099) q[1];
sx q[1];
rz(-0.0021541455) q[1];
rz(-pi) q[2];
rz(2.1625822) q[3];
sx q[3];
rz(-2.5755582) q[3];
sx q[3];
rz(2.4477959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4054823) q[2];
sx q[2];
rz(-2.5030899) q[2];
sx q[2];
rz(-0.68518266) q[2];
rz(-0.80125609) q[3];
sx q[3];
rz(-1.6359436) q[3];
sx q[3];
rz(-1.8906458) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41246978) q[0];
sx q[0];
rz(-2.6710489) q[0];
sx q[0];
rz(0.98010081) q[0];
rz(0.12067548) q[1];
sx q[1];
rz(-1.1680892) q[1];
sx q[1];
rz(1.4792222) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2187933) q[0];
sx q[0];
rz(-2.3258219) q[0];
sx q[0];
rz(-1.7936617) q[0];
rz(-pi) q[1];
rz(0.73320575) q[2];
sx q[2];
rz(-2.5016959) q[2];
sx q[2];
rz(-0.24240968) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4795717) q[1];
sx q[1];
rz(-1.5618854) q[1];
sx q[1];
rz(1.5655976) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9068662) q[3];
sx q[3];
rz(-1.0072636) q[3];
sx q[3];
rz(-0.8714217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4270619) q[2];
sx q[2];
rz(-1.4947299) q[2];
sx q[2];
rz(-0.18043268) q[2];
rz(0.42029542) q[3];
sx q[3];
rz(-0.97729483) q[3];
sx q[3];
rz(-0.54496566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9412398) q[0];
sx q[0];
rz(-1.7914597) q[0];
sx q[0];
rz(0.049276503) q[0];
rz(2.409528) q[1];
sx q[1];
rz(-2.2923636) q[1];
sx q[1];
rz(-1.3759618) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7017562) q[0];
sx q[0];
rz(-1.4539833) q[0];
sx q[0];
rz(-0.18420903) q[0];
rz(0.21993262) q[2];
sx q[2];
rz(-2.3637407) q[2];
sx q[2];
rz(-0.57524112) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0323309) q[1];
sx q[1];
rz(-1.9606661) q[1];
sx q[1];
rz(1.6435739) q[1];
x q[2];
rz(0.64582156) q[3];
sx q[3];
rz(-2.1649515) q[3];
sx q[3];
rz(0.37120968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0338318) q[2];
sx q[2];
rz(-1.7511448) q[2];
sx q[2];
rz(2.4793009) q[2];
rz(-1.3307339) q[3];
sx q[3];
rz(-0.8225421) q[3];
sx q[3];
rz(1.8432157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0419643) q[0];
sx q[0];
rz(-2.8350416) q[0];
sx q[0];
rz(-1.0076667) q[0];
rz(3.0362466) q[1];
sx q[1];
rz(-0.65709972) q[1];
sx q[1];
rz(-2.9109921) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4661575) q[0];
sx q[0];
rz(-0.49336067) q[0];
sx q[0];
rz(-2.3301279) q[0];
x q[1];
rz(-2.1449522) q[2];
sx q[2];
rz(-2.7189859) q[2];
sx q[2];
rz(-2.1017139) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9437792) q[1];
sx q[1];
rz(-0.18584419) q[1];
sx q[1];
rz(2.6617104) q[1];
x q[2];
rz(1.7218288) q[3];
sx q[3];
rz(-2.2103035) q[3];
sx q[3];
rz(-0.027627079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0838123) q[2];
sx q[2];
rz(-0.96692204) q[2];
sx q[2];
rz(-0.18873611) q[2];
rz(1.2356637) q[3];
sx q[3];
rz(-2.6344447) q[3];
sx q[3];
rz(-1.88131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96596232) q[0];
sx q[0];
rz(-1.2164793) q[0];
sx q[0];
rz(-2.8438582) q[0];
rz(-0.072602428) q[1];
sx q[1];
rz(-0.68521348) q[1];
sx q[1];
rz(2.4593478) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6202691) q[0];
sx q[0];
rz(-2.1730039) q[0];
sx q[0];
rz(-0.35767718) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14839006) q[2];
sx q[2];
rz(-2.1990657) q[2];
sx q[2];
rz(-0.65199696) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.025921019) q[1];
sx q[1];
rz(-2.762315) q[1];
sx q[1];
rz(1.161452) q[1];
rz(0.37118427) q[3];
sx q[3];
rz(-2.1371578) q[3];
sx q[3];
rz(-1.0584099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5785152) q[2];
sx q[2];
rz(-2.7882521) q[2];
sx q[2];
rz(-0.087470857) q[2];
rz(0.89720094) q[3];
sx q[3];
rz(-1.8950491) q[3];
sx q[3];
rz(2.7520666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0422269) q[0];
sx q[0];
rz(-2.0392188) q[0];
sx q[0];
rz(-1.1627831) q[0];
rz(-2.9258264) q[1];
sx q[1];
rz(-2.3240604) q[1];
sx q[1];
rz(-0.93528265) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1638086) q[0];
sx q[0];
rz(-0.80923128) q[0];
sx q[0];
rz(2.8092395) q[0];
x q[1];
rz(-1.8290524) q[2];
sx q[2];
rz(-1.4280945) q[2];
sx q[2];
rz(1.1772798) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.21970785) q[1];
sx q[1];
rz(-0.7898191) q[1];
sx q[1];
rz(0.11126693) q[1];
rz(-pi) q[2];
rz(2.8264826) q[3];
sx q[3];
rz(-1.2602046) q[3];
sx q[3];
rz(0.14652182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0995348) q[2];
sx q[2];
rz(-0.65379405) q[2];
sx q[2];
rz(0.81134861) q[2];
rz(-2.8387496) q[3];
sx q[3];
rz(-1.9767913) q[3];
sx q[3];
rz(2.5109049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.621156) q[0];
sx q[0];
rz(-2.8148837) q[0];
sx q[0];
rz(0.69123554) q[0];
rz(2.4825545) q[1];
sx q[1];
rz(-1.6233416) q[1];
sx q[1];
rz(2.5504327) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1189445) q[0];
sx q[0];
rz(-1.5592347) q[0];
sx q[0];
rz(-1.8537136) q[0];
rz(0.28556683) q[2];
sx q[2];
rz(-2.0519678) q[2];
sx q[2];
rz(2.8701631) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.65558479) q[1];
sx q[1];
rz(-0.84434768) q[1];
sx q[1];
rz(2.6996711) q[1];
rz(0.05686111) q[3];
sx q[3];
rz(-1.8627852) q[3];
sx q[3];
rz(2.9546471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9968694) q[2];
sx q[2];
rz(-2.150712) q[2];
sx q[2];
rz(1.6807751) q[2];
rz(0.75374323) q[3];
sx q[3];
rz(-1.4494579) q[3];
sx q[3];
rz(2.6678705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6853471) q[0];
sx q[0];
rz(-2.0401177) q[0];
sx q[0];
rz(-0.22612485) q[0];
rz(1.4489168) q[1];
sx q[1];
rz(-0.96335226) q[1];
sx q[1];
rz(2.1400145) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8346092) q[0];
sx q[0];
rz(-2.153647) q[0];
sx q[0];
rz(-1.7442763) q[0];
rz(1.1923157) q[2];
sx q[2];
rz(-2.2342302) q[2];
sx q[2];
rz(-2.0326234) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4521617) q[1];
sx q[1];
rz(-0.2161444) q[1];
sx q[1];
rz(2.6276905) q[1];
x q[2];
rz(-0.36205451) q[3];
sx q[3];
rz(-2.4241862) q[3];
sx q[3];
rz(-0.91594726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.61233026) q[2];
sx q[2];
rz(-0.77028209) q[2];
sx q[2];
rz(-0.15753499) q[2];
rz(-2.9869288) q[3];
sx q[3];
rz(-1.1636584) q[3];
sx q[3];
rz(-0.4304339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4100819) q[0];
sx q[0];
rz(-1.2393351) q[0];
sx q[0];
rz(2.4391158) q[0];
rz(-0.53600535) q[1];
sx q[1];
rz(-2.3119226) q[1];
sx q[1];
rz(-1.3943256) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047612543) q[0];
sx q[0];
rz(-1.5965921) q[0];
sx q[0];
rz(1.4071147) q[0];
rz(-pi) q[1];
rz(-0.093981103) q[2];
sx q[2];
rz(-1.6795484) q[2];
sx q[2];
rz(-0.55545872) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5750707) q[1];
sx q[1];
rz(-0.59986037) q[1];
sx q[1];
rz(-1.0529279) q[1];
rz(-1.5709747) q[3];
sx q[3];
rz(-1.9863673) q[3];
sx q[3];
rz(-0.10997575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1348306) q[2];
sx q[2];
rz(-0.68209058) q[2];
sx q[2];
rz(1.6949867) q[2];
rz(2.8805978) q[3];
sx q[3];
rz(-2.1311396) q[3];
sx q[3];
rz(-1.9058156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21988729) q[0];
sx q[0];
rz(-2.4867262) q[0];
sx q[0];
rz(2.126271) q[0];
rz(2.0657397) q[1];
sx q[1];
rz(-1.2747819) q[1];
sx q[1];
rz(-1.4631974) q[1];
rz(1.6816611) q[2];
sx q[2];
rz(-1.9456429) q[2];
sx q[2];
rz(2.6425998) q[2];
rz(-0.76293972) q[3];
sx q[3];
rz(-1.8992103) q[3];
sx q[3];
rz(-2.7873399) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
