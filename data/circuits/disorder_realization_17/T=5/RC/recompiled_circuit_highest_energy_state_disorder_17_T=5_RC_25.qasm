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
rz(-1.0626592) q[0];
sx q[0];
rz(-0.81544977) q[0];
sx q[0];
rz(0.71001473) q[0];
rz(2.8275936) q[1];
sx q[1];
rz(-2.2065838) q[1];
sx q[1];
rz(-1.8097872) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3478617) q[0];
sx q[0];
rz(-1.2958741) q[0];
sx q[0];
rz(1.3375651) q[0];
rz(-2.5153036) q[2];
sx q[2];
rz(-0.41538996) q[2];
sx q[2];
rz(2.9363869) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.59939042) q[1];
sx q[1];
rz(-2.1676461) q[1];
sx q[1];
rz(0.29920293) q[1];
rz(-1.0053357) q[3];
sx q[3];
rz(-0.71601235) q[3];
sx q[3];
rz(0.93040066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4189202) q[2];
sx q[2];
rz(-1.2166497) q[2];
sx q[2];
rz(-0.81532064) q[2];
rz(2.3319862) q[3];
sx q[3];
rz(-1.6428734) q[3];
sx q[3];
rz(2.0577551) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.078449) q[0];
sx q[0];
rz(-1.0493295) q[0];
sx q[0];
rz(0.11058841) q[0];
rz(2.1965006) q[1];
sx q[1];
rz(-0.4393622) q[1];
sx q[1];
rz(-1.5623215) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54540578) q[0];
sx q[0];
rz(-2.2431787) q[0];
sx q[0];
rz(0.59242146) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3198691) q[2];
sx q[2];
rz(-2.5603449) q[2];
sx q[2];
rz(1.7591214) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.64792127) q[1];
sx q[1];
rz(-2.646138) q[1];
sx q[1];
rz(0.22101553) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0840629) q[3];
sx q[3];
rz(-0.81765122) q[3];
sx q[3];
rz(0.14150894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.14757806) q[2];
sx q[2];
rz(-0.14573228) q[2];
sx q[2];
rz(-2.6914524) q[2];
rz(2.0077997) q[3];
sx q[3];
rz(-1.6885875) q[3];
sx q[3];
rz(-0.97761893) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56871539) q[0];
sx q[0];
rz(-0.94803888) q[0];
sx q[0];
rz(0.28133389) q[0];
rz(1.7549134) q[1];
sx q[1];
rz(-1.4298871) q[1];
sx q[1];
rz(-0.6001572) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2413797) q[0];
sx q[0];
rz(-0.29415392) q[0];
sx q[0];
rz(2.9808729) q[0];
rz(-0.05608989) q[2];
sx q[2];
rz(-0.96470913) q[2];
sx q[2];
rz(0.99585271) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0827209) q[1];
sx q[1];
rz(-1.2456254) q[1];
sx q[1];
rz(2.5356958) q[1];
x q[2];
rz(-0.26140499) q[3];
sx q[3];
rz(-2.5525744) q[3];
sx q[3];
rz(2.1847069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0448138) q[2];
sx q[2];
rz(-2.009095) q[2];
sx q[2];
rz(1.0463932) q[2];
rz(-2.7438296) q[3];
sx q[3];
rz(-0.64266509) q[3];
sx q[3];
rz(3.0090295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1716877) q[0];
sx q[0];
rz(-3.1145018) q[0];
sx q[0];
rz(1.0525674) q[0];
rz(-1.196208) q[1];
sx q[1];
rz(-1.2095249) q[1];
sx q[1];
rz(2.722091) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7561121) q[0];
sx q[0];
rz(-1.0051553) q[0];
sx q[0];
rz(1.2570981) q[0];
rz(-pi) q[1];
rz(-0.93650903) q[2];
sx q[2];
rz(-1.7542931) q[2];
sx q[2];
rz(-1.2681792) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0177815) q[1];
sx q[1];
rz(-2.3833125) q[1];
sx q[1];
rz(1.0067382) q[1];
x q[2];
rz(-2.64873) q[3];
sx q[3];
rz(-2.5774084) q[3];
sx q[3];
rz(-1.9436556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.18998751) q[2];
sx q[2];
rz(-2.0230468) q[2];
sx q[2];
rz(2.0513963) q[2];
rz(0.53127855) q[3];
sx q[3];
rz(-2.4254906) q[3];
sx q[3];
rz(1.5268911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7200318) q[0];
sx q[0];
rz(-1.5721385) q[0];
sx q[0];
rz(1.74362) q[0];
rz(-3.0086503) q[1];
sx q[1];
rz(-1.3016737) q[1];
sx q[1];
rz(-3.1403819) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.557363) q[0];
sx q[0];
rz(-1.5546636) q[0];
sx q[0];
rz(0.38059142) q[0];
rz(-pi) q[1];
rz(1.0693477) q[2];
sx q[2];
rz(-2.8184888) q[2];
sx q[2];
rz(0.94611514) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7442786) q[1];
sx q[1];
rz(-2.4220253) q[1];
sx q[1];
rz(-0.47081635) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4210971) q[3];
sx q[3];
rz(-1.5970917) q[3];
sx q[3];
rz(0.61248518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91263897) q[2];
sx q[2];
rz(-1.8043844) q[2];
sx q[2];
rz(-2.9310628) q[2];
rz(-2.1052965) q[3];
sx q[3];
rz(-2.3573037) q[3];
sx q[3];
rz(-1.9265296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1933111) q[0];
sx q[0];
rz(-1.223215) q[0];
sx q[0];
rz(2.7358828) q[0];
rz(-1.3544719) q[1];
sx q[1];
rz(-0.82258075) q[1];
sx q[1];
rz(-2.2192661) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.852823) q[0];
sx q[0];
rz(-1.6772224) q[0];
sx q[0];
rz(-1.2805416) q[0];
x q[1];
rz(3.0386557) q[2];
sx q[2];
rz(-0.90667001) q[2];
sx q[2];
rz(0.31598202) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3971484) q[1];
sx q[1];
rz(-1.4566629) q[1];
sx q[1];
rz(0.95790095) q[1];
x q[2];
rz(-2.5821463) q[3];
sx q[3];
rz(-1.2370951) q[3];
sx q[3];
rz(-2.4041686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7102082) q[2];
sx q[2];
rz(-1.9269383) q[2];
sx q[2];
rz(2.7396438) q[2];
rz(-1.8534144) q[3];
sx q[3];
rz(-0.86789075) q[3];
sx q[3];
rz(1.2791546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5227018) q[0];
sx q[0];
rz(-1.1033449) q[0];
sx q[0];
rz(0.064362137) q[0];
rz(0.83802682) q[1];
sx q[1];
rz(-0.82632724) q[1];
sx q[1];
rz(1.3305957) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2272328) q[0];
sx q[0];
rz(-1.6961251) q[0];
sx q[0];
rz(-2.5663239) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0688017) q[2];
sx q[2];
rz(-1.0235707) q[2];
sx q[2];
rz(-0.71820532) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.76112642) q[1];
sx q[1];
rz(-2.5674324) q[1];
sx q[1];
rz(2.7284184) q[1];
x q[2];
rz(1.3279541) q[3];
sx q[3];
rz(-2.5153219) q[3];
sx q[3];
rz(2.0988718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7360721) q[2];
sx q[2];
rz(-1.5084927) q[2];
sx q[2];
rz(2.4580477) q[2];
rz(0.73379597) q[3];
sx q[3];
rz(-1.4132696) q[3];
sx q[3];
rz(-2.2939513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7965294) q[0];
sx q[0];
rz(-2.0579484) q[0];
sx q[0];
rz(-1.9684568) q[0];
rz(0.37551156) q[1];
sx q[1];
rz(-2.651732) q[1];
sx q[1];
rz(-2.1048022) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11060729) q[0];
sx q[0];
rz(-1.5944423) q[0];
sx q[0];
rz(-1.5899237) q[0];
x q[1];
rz(2.1027742) q[2];
sx q[2];
rz(-0.3178645) q[2];
sx q[2];
rz(-1.7828072) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.97809659) q[1];
sx q[1];
rz(-1.4304203) q[1];
sx q[1];
rz(2.3771493) q[1];
x q[2];
rz(2.6508207) q[3];
sx q[3];
rz(-2.0476855) q[3];
sx q[3];
rz(-0.27974883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5658687) q[2];
sx q[2];
rz(-0.91369358) q[2];
sx q[2];
rz(-2.4110528) q[2];
rz(-0.20914397) q[3];
sx q[3];
rz(-2.6181965) q[3];
sx q[3];
rz(1.2701579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4850979) q[0];
sx q[0];
rz(-0.42335835) q[0];
sx q[0];
rz(-3.094161) q[0];
rz(2.9196396) q[1];
sx q[1];
rz(-2.0675979) q[1];
sx q[1];
rz(2.2304631) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.297356) q[0];
sx q[0];
rz(-0.043023303) q[0];
sx q[0];
rz(0.6776612) q[0];
x q[1];
rz(1.59052) q[2];
sx q[2];
rz(-2.4998186) q[2];
sx q[2];
rz(0.80427792) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4621689) q[1];
sx q[1];
rz(-2.7738681) q[1];
sx q[1];
rz(3.1064139) q[1];
rz(-pi) q[2];
rz(0.20424517) q[3];
sx q[3];
rz(-2.1276655) q[3];
sx q[3];
rz(2.000252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.86965108) q[2];
sx q[2];
rz(-0.38001529) q[2];
sx q[2];
rz(-0.16858777) q[2];
rz(-2.9224959) q[3];
sx q[3];
rz(-2.2897661) q[3];
sx q[3];
rz(1.3271837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1472226) q[0];
sx q[0];
rz(-2.8012025) q[0];
sx q[0];
rz(-0.45743531) q[0];
rz(-1.9153197) q[1];
sx q[1];
rz(-2.8044082) q[1];
sx q[1];
rz(-1.3630684) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6493312) q[0];
sx q[0];
rz(-2.1773585) q[0];
sx q[0];
rz(-1.6442293) q[0];
rz(-pi) q[1];
rz(2.5117433) q[2];
sx q[2];
rz(-1.7694663) q[2];
sx q[2];
rz(0.82306403) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5879092) q[1];
sx q[1];
rz(-1.739555) q[1];
sx q[1];
rz(-0.38270271) q[1];
rz(0.8882723) q[3];
sx q[3];
rz(-1.5049045) q[3];
sx q[3];
rz(0.94277387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3739796) q[2];
sx q[2];
rz(-2.5869936) q[2];
sx q[2];
rz(2.6895831) q[2];
rz(-2.7376145) q[3];
sx q[3];
rz(-1.2510866) q[3];
sx q[3];
rz(-1.1261136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44611888) q[0];
sx q[0];
rz(-1.5344545) q[0];
sx q[0];
rz(-2.9547966) q[0];
rz(-0.52234621) q[1];
sx q[1];
rz(-0.53032395) q[1];
sx q[1];
rz(1.2293336) q[1];
rz(-0.73595388) q[2];
sx q[2];
rz(-2.7685168) q[2];
sx q[2];
rz(-1.4442486) q[2];
rz(-2.8780739) q[3];
sx q[3];
rz(-1.4265887) q[3];
sx q[3];
rz(-1.9909457) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
