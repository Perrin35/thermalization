OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2394387) q[0];
sx q[0];
rz(-1.8128938) q[0];
sx q[0];
rz(-0.21188307) q[0];
rz(-1.4172685) q[1];
sx q[1];
rz(-2.6091726) q[1];
sx q[1];
rz(0.37791696) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48403063) q[0];
sx q[0];
rz(-3.1320509) q[0];
sx q[0];
rz(1.5242759) q[0];
rz(2.7013999) q[2];
sx q[2];
rz(-0.72840103) q[2];
sx q[2];
rz(0.0022526646) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.35740543) q[1];
sx q[1];
rz(-2.2012246) q[1];
sx q[1];
rz(1.9131843) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.7805232) q[3];
sx q[3];
rz(-0.98376432) q[3];
sx q[3];
rz(-0.39654708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2892896) q[2];
sx q[2];
rz(-2.0630344) q[2];
sx q[2];
rz(-3.0618073) q[2];
rz(-2.1763109) q[3];
sx q[3];
rz(-1.691317) q[3];
sx q[3];
rz(0.7307581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6642283) q[0];
sx q[0];
rz(-2.1381162) q[0];
sx q[0];
rz(0.088951237) q[0];
rz(1.7695919) q[1];
sx q[1];
rz(-1.7141432) q[1];
sx q[1];
rz(2.037183) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0029582214) q[0];
sx q[0];
rz(-1.3413652) q[0];
sx q[0];
rz(-2.5087439) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78271336) q[2];
sx q[2];
rz(-1.209895) q[2];
sx q[2];
rz(-0.95907839) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4648393) q[1];
sx q[1];
rz(-2.0130806) q[1];
sx q[1];
rz(-0.16525903) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9761451) q[3];
sx q[3];
rz(-2.5962127) q[3];
sx q[3];
rz(-0.65217962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8770807) q[2];
sx q[2];
rz(-2.0161714) q[2];
sx q[2];
rz(1.4667) q[2];
rz(-3.1125715) q[3];
sx q[3];
rz(-2.1252188) q[3];
sx q[3];
rz(-2.4513054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2324209) q[0];
sx q[0];
rz(-3.0747774) q[0];
sx q[0];
rz(-2.8636041) q[0];
rz(1.7104644) q[1];
sx q[1];
rz(-2.2093096) q[1];
sx q[1];
rz(2.7412282) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5059197) q[0];
sx q[0];
rz(-2.189496) q[0];
sx q[0];
rz(2.1257504) q[0];
rz(0.68636559) q[2];
sx q[2];
rz(-1.7698506) q[2];
sx q[2];
rz(0.65716568) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0644933) q[1];
sx q[1];
rz(-2.0410186) q[1];
sx q[1];
rz(1.749586) q[1];
rz(-1.2648029) q[3];
sx q[3];
rz(-1.8974432) q[3];
sx q[3];
rz(-1.6040925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4572738) q[2];
sx q[2];
rz(-1.6088586) q[2];
sx q[2];
rz(2.686783) q[2];
rz(-1.2416035) q[3];
sx q[3];
rz(-2.1865215) q[3];
sx q[3];
rz(0.72030592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1860745) q[0];
sx q[0];
rz(-1.3294514) q[0];
sx q[0];
rz(-0.00051001471) q[0];
rz(-2.5406802) q[1];
sx q[1];
rz(-2.3020703) q[1];
sx q[1];
rz(-2.9972163) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52989456) q[0];
sx q[0];
rz(-1.7774044) q[0];
sx q[0];
rz(2.1312461) q[0];
rz(-1.7208485) q[2];
sx q[2];
rz(-1.320082) q[2];
sx q[2];
rz(-1.5992407) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.51203242) q[1];
sx q[1];
rz(-1.343439) q[1];
sx q[1];
rz(-2.4870706) q[1];
rz(-pi) q[2];
rz(0.17474971) q[3];
sx q[3];
rz(-0.55395444) q[3];
sx q[3];
rz(0.21077158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1986177) q[2];
sx q[2];
rz(-1.4451507) q[2];
sx q[2];
rz(0.96251881) q[2];
rz(1.5396384) q[3];
sx q[3];
rz(-1.4055777) q[3];
sx q[3];
rz(-0.34077728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9050423) q[0];
sx q[0];
rz(-1.4555229) q[0];
sx q[0];
rz(1.1849674) q[0];
rz(-0.22625893) q[1];
sx q[1];
rz(-0.87892756) q[1];
sx q[1];
rz(2.102898) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13030355) q[0];
sx q[0];
rz(-1.652392) q[0];
sx q[0];
rz(-2.5528583) q[0];
rz(-pi) q[1];
rz(2.2009497) q[2];
sx q[2];
rz(-2.2046208) q[2];
sx q[2];
rz(1.188082) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2782368) q[1];
sx q[1];
rz(-1.0627295) q[1];
sx q[1];
rz(2.3418576) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4323147) q[3];
sx q[3];
rz(-2.6977728) q[3];
sx q[3];
rz(0.5139155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.431939) q[2];
sx q[2];
rz(-0.69810549) q[2];
sx q[2];
rz(-2.8254438) q[2];
rz(1.4656674) q[3];
sx q[3];
rz(-1.1301872) q[3];
sx q[3];
rz(-1.7620311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6370711) q[0];
sx q[0];
rz(-2.2383454) q[0];
sx q[0];
rz(2.0380518) q[0];
rz(-2.6612813) q[1];
sx q[1];
rz(-2.4791398) q[1];
sx q[1];
rz(-0.59741098) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53086583) q[0];
sx q[0];
rz(-1.382575) q[0];
sx q[0];
rz(-1.782531) q[0];
rz(0.74540794) q[2];
sx q[2];
rz(-1.3920203) q[2];
sx q[2];
rz(2.488236) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0805507) q[1];
sx q[1];
rz(-1.7361199) q[1];
sx q[1];
rz(-0.25456659) q[1];
rz(-1.3378998) q[3];
sx q[3];
rz(-0.38023708) q[3];
sx q[3];
rz(-1.5986795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.20208134) q[2];
sx q[2];
rz(-1.6263522) q[2];
sx q[2];
rz(-2.0651979) q[2];
rz(1.1427897) q[3];
sx q[3];
rz(-2.3484774) q[3];
sx q[3];
rz(-2.6514261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48653212) q[0];
sx q[0];
rz(-1.9831816) q[0];
sx q[0];
rz(3.1031188) q[0];
rz(-0.064727457) q[1];
sx q[1];
rz(-1.3836626) q[1];
sx q[1];
rz(0.23385349) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9486987) q[0];
sx q[0];
rz(-1.3854376) q[0];
sx q[0];
rz(2.9167487) q[0];
rz(-pi) q[1];
rz(2.269417) q[2];
sx q[2];
rz(-2.6769014) q[2];
sx q[2];
rz(1.8605491) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.13713) q[1];
sx q[1];
rz(-0.65441583) q[1];
sx q[1];
rz(0.14738247) q[1];
rz(-pi) q[2];
rz(3.0117118) q[3];
sx q[3];
rz(-1.5135153) q[3];
sx q[3];
rz(-0.86806017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6528299) q[2];
sx q[2];
rz(-1.5583928) q[2];
sx q[2];
rz(-0.59824198) q[2];
rz(-3.0277142) q[3];
sx q[3];
rz(-1.7555534) q[3];
sx q[3];
rz(0.84754506) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7365731) q[0];
sx q[0];
rz(-2.2990655) q[0];
sx q[0];
rz(2.8952428) q[0];
rz(1.8440638) q[1];
sx q[1];
rz(-1.1187226) q[1];
sx q[1];
rz(-2.6447703) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4637359) q[0];
sx q[0];
rz(-1.1564768) q[0];
sx q[0];
rz(0.68092771) q[0];
rz(2.970817) q[2];
sx q[2];
rz(-2.429212) q[2];
sx q[2];
rz(-2.8757489) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0973952) q[1];
sx q[1];
rz(-2.5977511) q[1];
sx q[1];
rz(1.103765) q[1];
rz(-pi) q[2];
rz(1.0786177) q[3];
sx q[3];
rz(-1.4340622) q[3];
sx q[3];
rz(-2.0355952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7184427) q[2];
sx q[2];
rz(-2.6079874) q[2];
sx q[2];
rz(-1.9971087) q[2];
rz(-1.9874969) q[3];
sx q[3];
rz(-1.4242947) q[3];
sx q[3];
rz(0.19788876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4474739) q[0];
sx q[0];
rz(-0.23451528) q[0];
sx q[0];
rz(-3.0754454) q[0];
rz(-1.9937531) q[1];
sx q[1];
rz(-1.3775974) q[1];
sx q[1];
rz(2.537421) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2691374) q[0];
sx q[0];
rz(-0.34664422) q[0];
sx q[0];
rz(-2.4983062) q[0];
rz(-0.2653052) q[2];
sx q[2];
rz(-0.54232208) q[2];
sx q[2];
rz(0.10077439) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.79163247) q[1];
sx q[1];
rz(-2.5436328) q[1];
sx q[1];
rz(2.7954742) q[1];
rz(0.87852134) q[3];
sx q[3];
rz(-2.1067348) q[3];
sx q[3];
rz(1.5978447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8784647) q[2];
sx q[2];
rz(-0.85004127) q[2];
sx q[2];
rz(2.7086332) q[2];
rz(1.2290907) q[3];
sx q[3];
rz(-1.9357598) q[3];
sx q[3];
rz(1.3273201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1935254) q[0];
sx q[0];
rz(-2.0663517) q[0];
sx q[0];
rz(-1.2731592) q[0];
rz(0.46514568) q[1];
sx q[1];
rz(-1.3659313) q[1];
sx q[1];
rz(-2.926631) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0996159) q[0];
sx q[0];
rz(-0.93659217) q[0];
sx q[0];
rz(-1.8862104) q[0];
rz(-pi) q[1];
rz(-0.6647756) q[2];
sx q[2];
rz(-0.6419581) q[2];
sx q[2];
rz(2.4017815) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.82483208) q[1];
sx q[1];
rz(-2.2344898) q[1];
sx q[1];
rz(0.65726991) q[1];
rz(-pi) q[2];
rz(-0.56193476) q[3];
sx q[3];
rz(-1.7896381) q[3];
sx q[3];
rz(-2.9666025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.8421858) q[2];
sx q[2];
rz(-0.81373787) q[2];
sx q[2];
rz(-2.1827533) q[2];
rz(1.9793319) q[3];
sx q[3];
rz(-1.6543038) q[3];
sx q[3];
rz(1.4452665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5398298) q[0];
sx q[0];
rz(-1.6015263) q[0];
sx q[0];
rz(1.482561) q[0];
rz(2.4330347) q[1];
sx q[1];
rz(-0.19150145) q[1];
sx q[1];
rz(2.3932744) q[1];
rz(-0.49025771) q[2];
sx q[2];
rz(-0.62050642) q[2];
sx q[2];
rz(-1.4668589) q[2];
rz(1.7583354) q[3];
sx q[3];
rz(-2.2839727) q[3];
sx q[3];
rz(-1.4061389) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
