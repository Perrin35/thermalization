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
rz(3.6291549) q[0];
sx q[0];
rz(3.5222375) q[0];
sx q[0];
rz(9.335523) q[0];
rz(-3.6042876) q[1];
sx q[1];
rz(2.8609639) q[1];
sx q[1];
rz(14.345471) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93264047) q[0];
sx q[0];
rz(-2.7679043) q[0];
sx q[0];
rz(-3.0840918) q[0];
x q[1];
rz(1.5342185) q[2];
sx q[2];
rz(-1.4223411) q[2];
sx q[2];
rz(-1.0561933) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.57682788) q[1];
sx q[1];
rz(-2.514808) q[1];
sx q[1];
rz(0.38615055) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6854679) q[3];
sx q[3];
rz(-1.4906192) q[3];
sx q[3];
rz(-2.2661569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6224299) q[2];
sx q[2];
rz(-0.76202718) q[2];
sx q[2];
rz(-2.2117174) q[2];
rz(0.30396384) q[3];
sx q[3];
rz(-1.8094939) q[3];
sx q[3];
rz(-2.9728594) q[3];
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
rz(-pi) q[0];
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
rz(0.25794432) q[0];
sx q[0];
rz(-1.5644263) q[0];
sx q[0];
rz(1.2516578) q[0];
rz(-0.21580639) q[1];
sx q[1];
rz(-1.0390751) q[1];
sx q[1];
rz(3.0024517) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9342095) q[0];
sx q[0];
rz(-1.73044) q[0];
sx q[0];
rz(-0.90682323) q[0];
rz(-0.49824841) q[2];
sx q[2];
rz(-1.3300657) q[2];
sx q[2];
rz(1.6290851) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5522668) q[1];
sx q[1];
rz(-1.3331474) q[1];
sx q[1];
rz(0.71961211) q[1];
rz(-pi) q[2];
rz(0.82461951) q[3];
sx q[3];
rz(-2.1420205) q[3];
sx q[3];
rz(-2.415641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4591878) q[2];
sx q[2];
rz(-2.6725957) q[2];
sx q[2];
rz(-1.0914618) q[2];
rz(-1.5952716) q[3];
sx q[3];
rz(-1.073444) q[3];
sx q[3];
rz(0.44927621) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0794373) q[0];
sx q[0];
rz(-1.2261483) q[0];
sx q[0];
rz(-0.011818258) q[0];
rz(-0.91105175) q[1];
sx q[1];
rz(-1.3900737) q[1];
sx q[1];
rz(-1.3284838) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8810711) q[0];
sx q[0];
rz(-2.0876472) q[0];
sx q[0];
rz(-1.757574) q[0];
rz(-pi) q[1];
rz(-2.6027868) q[2];
sx q[2];
rz(-2.9201815) q[2];
sx q[2];
rz(3.0543229) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3551533) q[1];
sx q[1];
rz(-2.3117073) q[1];
sx q[1];
rz(2.607858) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4592917) q[3];
sx q[3];
rz(-1.6435677) q[3];
sx q[3];
rz(2.7804216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7517684) q[2];
sx q[2];
rz(-1.5247034) q[2];
sx q[2];
rz(-0.41979182) q[2];
rz(-1.3587562) q[3];
sx q[3];
rz(-0.32521453) q[3];
sx q[3];
rz(1.1903919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(-1.457394) q[0];
sx q[0];
rz(-2.0036819) q[0];
sx q[0];
rz(-2.6275291) q[0];
rz(-0.38814107) q[1];
sx q[1];
rz(-0.61182794) q[1];
sx q[1];
rz(1.911389) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82321754) q[0];
sx q[0];
rz(-1.3127232) q[0];
sx q[0];
rz(2.8283872) q[0];
rz(-0.72429652) q[2];
sx q[2];
rz(-1.8593374) q[2];
sx q[2];
rz(-2.9674825) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9446805) q[1];
sx q[1];
rz(-0.80183235) q[1];
sx q[1];
rz(-0.51735984) q[1];
x q[2];
rz(-0.28546412) q[3];
sx q[3];
rz(-2.450911) q[3];
sx q[3];
rz(2.635298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6396883) q[2];
sx q[2];
rz(-1.014726) q[2];
sx q[2];
rz(-0.91628966) q[2];
rz(-2.6472951) q[3];
sx q[3];
rz(-1.1980134) q[3];
sx q[3];
rz(2.3673207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2275527) q[0];
sx q[0];
rz(-0.95714772) q[0];
sx q[0];
rz(-0.48602948) q[0];
rz(0.60274974) q[1];
sx q[1];
rz(-1.175368) q[1];
sx q[1];
rz(1.1154307) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91624068) q[0];
sx q[0];
rz(-1.1778206) q[0];
sx q[0];
rz(2.4797702) q[0];
x q[1];
rz(2.7390295) q[2];
sx q[2];
rz(-0.83588119) q[2];
sx q[2];
rz(-2.1083567) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.095905) q[1];
sx q[1];
rz(-2.3153911) q[1];
sx q[1];
rz(0.5788486) q[1];
rz(-2.0492184) q[3];
sx q[3];
rz(-0.53077519) q[3];
sx q[3];
rz(-1.9984744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9976161) q[2];
sx q[2];
rz(-2.82085) q[2];
sx q[2];
rz(2.0799267) q[2];
rz(1.5205421) q[3];
sx q[3];
rz(-1.942626) q[3];
sx q[3];
rz(-2.2475713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(3.0735556) q[0];
sx q[0];
rz(-3.1017922) q[0];
sx q[0];
rz(-1.918248) q[0];
rz(-2.6262737) q[1];
sx q[1];
rz(-0.66910187) q[1];
sx q[1];
rz(1.7759148) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0037352) q[0];
sx q[0];
rz(-0.23085871) q[0];
sx q[0];
rz(-1.0539216) q[0];
x q[1];
rz(3.0480644) q[2];
sx q[2];
rz(-1.3892738) q[2];
sx q[2];
rz(-2.7101171) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.15348831) q[1];
sx q[1];
rz(-2.0447013) q[1];
sx q[1];
rz(-1.7405645) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29748282) q[3];
sx q[3];
rz(-1.6053891) q[3];
sx q[3];
rz(-2.5643947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0627275) q[2];
sx q[2];
rz(-1.8026423) q[2];
sx q[2];
rz(-1.3272746) q[2];
rz(-1.7257388) q[3];
sx q[3];
rz(-3.0684107) q[3];
sx q[3];
rz(-2.2590526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0075204) q[0];
sx q[0];
rz(-0.53711689) q[0];
sx q[0];
rz(2.9448331) q[0];
rz(-0.21601954) q[1];
sx q[1];
rz(-1.0088423) q[1];
sx q[1];
rz(-0.78741995) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24889937) q[0];
sx q[0];
rz(-1.8149879) q[0];
sx q[0];
rz(-2.9466183) q[0];
x q[1];
rz(0.51047275) q[2];
sx q[2];
rz(-1.831448) q[2];
sx q[2];
rz(2.5939121) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.43147165) q[1];
sx q[1];
rz(-0.97968757) q[1];
sx q[1];
rz(1.6256871) q[1];
x q[2];
rz(-2.092497) q[3];
sx q[3];
rz(-0.36715436) q[3];
sx q[3];
rz(3.1279813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6009723) q[2];
sx q[2];
rz(-1.5121907) q[2];
sx q[2];
rz(-2.7579894) q[2];
rz(-2.2331734) q[3];
sx q[3];
rz(-0.81276613) q[3];
sx q[3];
rz(2.9878591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8613794) q[0];
sx q[0];
rz(-2.6160243) q[0];
sx q[0];
rz(0.63661611) q[0];
rz(-0.55066291) q[1];
sx q[1];
rz(-1.5148342) q[1];
sx q[1];
rz(0.47264636) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022332683) q[0];
sx q[0];
rz(-0.68281931) q[0];
sx q[0];
rz(-0.14108087) q[0];
x q[1];
rz(1.6238975) q[2];
sx q[2];
rz(-1.6638954) q[2];
sx q[2];
rz(2.0123717) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0090985) q[1];
sx q[1];
rz(-0.14202296) q[1];
sx q[1];
rz(0.36722398) q[1];
x q[2];
rz(1.8697463) q[3];
sx q[3];
rz(-0.96072324) q[3];
sx q[3];
rz(-1.7536193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7134646) q[2];
sx q[2];
rz(-1.8992004) q[2];
sx q[2];
rz(0.098527519) q[2];
rz(-0.030700961) q[3];
sx q[3];
rz(-1.5714329) q[3];
sx q[3];
rz(2.9876685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6105662) q[0];
sx q[0];
rz(-2.1838146) q[0];
sx q[0];
rz(-0.73262334) q[0];
rz(-3.1351807) q[1];
sx q[1];
rz(-2.1156204) q[1];
sx q[1];
rz(1.7195255) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61056449) q[0];
sx q[0];
rz(-2.8262666) q[0];
sx q[0];
rz(0.91489961) q[0];
x q[1];
rz(-1.9344011) q[2];
sx q[2];
rz(-1.0997314) q[2];
sx q[2];
rz(-2.1944012) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.040557794) q[1];
sx q[1];
rz(-0.15753105) q[1];
sx q[1];
rz(-2.9228802) q[1];
x q[2];
rz(-0.51958618) q[3];
sx q[3];
rz(-1.6359463) q[3];
sx q[3];
rz(1.6217854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0942568) q[2];
sx q[2];
rz(-1.3553268) q[2];
sx q[2];
rz(-2.6195841) q[2];
rz(-3.1212854) q[3];
sx q[3];
rz(-2.4419407) q[3];
sx q[3];
rz(0.29683963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8727528) q[0];
sx q[0];
rz(-1.5248542) q[0];
sx q[0];
rz(-2.010345) q[0];
rz(-2.3616683) q[1];
sx q[1];
rz(-1.9386539) q[1];
sx q[1];
rz(-2.6663229) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4085379) q[0];
sx q[0];
rz(-2.2768639) q[0];
sx q[0];
rz(0.48581328) q[0];
rz(2.2384727) q[2];
sx q[2];
rz(-1.3851056) q[2];
sx q[2];
rz(-2.8201495) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1058029) q[1];
sx q[1];
rz(-2.8052108) q[1];
sx q[1];
rz(0.85806429) q[1];
rz(-pi) q[2];
rz(-0.69644955) q[3];
sx q[3];
rz(-1.4587194) q[3];
sx q[3];
rz(-0.96986412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7705226) q[2];
sx q[2];
rz(-1.5757685) q[2];
sx q[2];
rz(-1.9353665) q[2];
rz(1.3953588) q[3];
sx q[3];
rz(-0.54308707) q[3];
sx q[3];
rz(3.0622845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64410011) q[0];
sx q[0];
rz(-2.946749) q[0];
sx q[0];
rz(-0.83644833) q[0];
rz(-2.1438228) q[1];
sx q[1];
rz(-1.0918959) q[1];
sx q[1];
rz(0.087654884) q[1];
rz(1.8733415) q[2];
sx q[2];
rz(-0.22845636) q[2];
sx q[2];
rz(-0.34194389) q[2];
rz(-1.6818123) q[3];
sx q[3];
rz(-0.89499499) q[3];
sx q[3];
rz(0.53193308) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
