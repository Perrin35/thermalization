OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0290282) q[0];
sx q[0];
rz(-1.5052786) q[0];
sx q[0];
rz(2.3583052) q[0];
rz(-2.0593491) q[1];
sx q[1];
rz(-1.6545656) q[1];
sx q[1];
rz(2.3496871) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5675674) q[0];
sx q[0];
rz(-1.3114531) q[0];
sx q[0];
rz(-1.7631084) q[0];
rz(-1.7246805) q[2];
sx q[2];
rz(-1.6072011) q[2];
sx q[2];
rz(-2.3652181) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.82552443) q[1];
sx q[1];
rz(-1.5393942) q[1];
sx q[1];
rz(-2.5981748) q[1];
rz(0.79444076) q[3];
sx q[3];
rz(-2.0047965) q[3];
sx q[3];
rz(1.7931012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.59114328) q[2];
sx q[2];
rz(-0.14269665) q[2];
sx q[2];
rz(0.099451065) q[2];
rz(0.24672306) q[3];
sx q[3];
rz(-1.6171425) q[3];
sx q[3];
rz(2.7439086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1219015) q[0];
sx q[0];
rz(-2.1653403) q[0];
sx q[0];
rz(-0.9285399) q[0];
rz(1.4506725) q[1];
sx q[1];
rz(-1.848315) q[1];
sx q[1];
rz(0.64847747) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6715393) q[0];
sx q[0];
rz(-2.3190365) q[0];
sx q[0];
rz(2.3683192) q[0];
rz(1.2546059) q[2];
sx q[2];
rz(-1.0266227) q[2];
sx q[2];
rz(2.3426825) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.14997031) q[1];
sx q[1];
rz(-0.14279616) q[1];
sx q[1];
rz(0.10985978) q[1];
rz(-1.2821372) q[3];
sx q[3];
rz(-1.4573754) q[3];
sx q[3];
rz(0.78383776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4497946) q[2];
sx q[2];
rz(-2.3851676) q[2];
sx q[2];
rz(2.2890384) q[2];
rz(-2.2954588) q[3];
sx q[3];
rz(-1.5087912) q[3];
sx q[3];
rz(-1.030863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6344675) q[0];
sx q[0];
rz(-1.9159303) q[0];
sx q[0];
rz(1.0852098) q[0];
rz(2.197544) q[1];
sx q[1];
rz(-1.4986821) q[1];
sx q[1];
rz(1.5012213) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47848343) q[0];
sx q[0];
rz(-1.5318277) q[0];
sx q[0];
rz(-2.9443254) q[0];
rz(-pi) q[1];
rz(-2.2527462) q[2];
sx q[2];
rz(-0.29817981) q[2];
sx q[2];
rz(-2.2928638) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2861917) q[1];
sx q[1];
rz(-2.5603189) q[1];
sx q[1];
rz(1.3894677) q[1];
rz(-0.94475475) q[3];
sx q[3];
rz(-0.26492719) q[3];
sx q[3];
rz(-2.0719599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7368855) q[2];
sx q[2];
rz(-0.53386226) q[2];
sx q[2];
rz(0.85477465) q[2];
rz(-2.2331623) q[3];
sx q[3];
rz(-1.5646489) q[3];
sx q[3];
rz(-2.0250208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.214355) q[0];
sx q[0];
rz(-2.7404009) q[0];
sx q[0];
rz(1.9450564) q[0];
rz(-1.6664956) q[1];
sx q[1];
rz(-2.7207082) q[1];
sx q[1];
rz(-1.9570785) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9036983) q[0];
sx q[0];
rz(-1.333295) q[0];
sx q[0];
rz(-2.4303275) q[0];
rz(-pi) q[1];
x q[1];
rz(2.938561) q[2];
sx q[2];
rz(-0.76329279) q[2];
sx q[2];
rz(-0.53084669) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.056442913) q[1];
sx q[1];
rz(-0.55906017) q[1];
sx q[1];
rz(2.0075625) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7538025) q[3];
sx q[3];
rz(-1.840045) q[3];
sx q[3];
rz(1.315801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2148332) q[2];
sx q[2];
rz(-0.49761179) q[2];
sx q[2];
rz(1.8801749) q[2];
rz(2.7091806) q[3];
sx q[3];
rz(-1.132248) q[3];
sx q[3];
rz(2.6692218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.072902) q[0];
sx q[0];
rz(-1.9111159) q[0];
sx q[0];
rz(3.1121837) q[0];
rz(-0.75621653) q[1];
sx q[1];
rz(-0.58719802) q[1];
sx q[1];
rz(-0.75278935) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8653976) q[0];
sx q[0];
rz(-1.5711693) q[0];
sx q[0];
rz(1.5764109) q[0];
x q[1];
rz(1.8389614) q[2];
sx q[2];
rz(-1.6577621) q[2];
sx q[2];
rz(0.11448569) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.50548762) q[1];
sx q[1];
rz(-1.6754181) q[1];
sx q[1];
rz(-0.57144798) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1208956) q[3];
sx q[3];
rz(-2.0664762) q[3];
sx q[3];
rz(-0.97199398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0011562) q[2];
sx q[2];
rz(-1.493528) q[2];
sx q[2];
rz(2.746554) q[2];
rz(0.47518528) q[3];
sx q[3];
rz(-1.7893712) q[3];
sx q[3];
rz(2.7648259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7659371) q[0];
sx q[0];
rz(-0.7203311) q[0];
sx q[0];
rz(0.96250594) q[0];
rz(-0.51482254) q[1];
sx q[1];
rz(-1.5341026) q[1];
sx q[1];
rz(2.8033676) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048934919) q[0];
sx q[0];
rz(-0.83507628) q[0];
sx q[0];
rz(0.2661163) q[0];
x q[1];
rz(-1.3231311) q[2];
sx q[2];
rz(-1.3405352) q[2];
sx q[2];
rz(2.3375653) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.64778642) q[1];
sx q[1];
rz(-1.5739417) q[1];
sx q[1];
rz(3.1415658) q[1];
x q[2];
rz(1.4184444) q[3];
sx q[3];
rz(-0.42644106) q[3];
sx q[3];
rz(-1.3063198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.50011355) q[2];
sx q[2];
rz(-1.1621472) q[2];
sx q[2];
rz(2.5872453) q[2];
rz(0.85136271) q[3];
sx q[3];
rz(-2.7832289) q[3];
sx q[3];
rz(0.62801492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8165269) q[0];
sx q[0];
rz(-2.5287703) q[0];
sx q[0];
rz(-2.391173) q[0];
rz(2.5665307) q[1];
sx q[1];
rz(-1.776418) q[1];
sx q[1];
rz(-0.65779984) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3391425) q[0];
sx q[0];
rz(-0.93451148) q[0];
sx q[0];
rz(-2.0872981) q[0];
rz(3.0900938) q[2];
sx q[2];
rz(-2.4337075) q[2];
sx q[2];
rz(-2.5619363) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.15193394) q[1];
sx q[1];
rz(-0.2865782) q[1];
sx q[1];
rz(-1.5119988) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9255387) q[3];
sx q[3];
rz(-1.1974058) q[3];
sx q[3];
rz(1.7301529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.49360069) q[2];
sx q[2];
rz(-2.1330264) q[2];
sx q[2];
rz(-1.4823401) q[2];
rz(0.12065398) q[3];
sx q[3];
rz(-1.7574661) q[3];
sx q[3];
rz(0.83546662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3487314) q[0];
sx q[0];
rz(-2.272235) q[0];
sx q[0];
rz(-0.29801512) q[0];
rz(2.1030078) q[1];
sx q[1];
rz(-0.54439259) q[1];
sx q[1];
rz(0.15377741) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34798056) q[0];
sx q[0];
rz(-1.4958994) q[0];
sx q[0];
rz(2.8872284) q[0];
rz(-0.80259097) q[2];
sx q[2];
rz(-1.0804515) q[2];
sx q[2];
rz(-2.5732694) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.41339918) q[1];
sx q[1];
rz(-1.4811885) q[1];
sx q[1];
rz(2.7941568) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12184398) q[3];
sx q[3];
rz(-1.9949556) q[3];
sx q[3];
rz(0.76790732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.095857233) q[2];
sx q[2];
rz(-0.46892527) q[2];
sx q[2];
rz(-1.1886965) q[2];
rz(-1.3304322) q[3];
sx q[3];
rz(-1.9537787) q[3];
sx q[3];
rz(2.4082898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34340149) q[0];
sx q[0];
rz(-0.60438406) q[0];
sx q[0];
rz(-1.6424302) q[0];
rz(2.5921953) q[1];
sx q[1];
rz(-1.5769985) q[1];
sx q[1];
rz(-2.8909491) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9662387) q[0];
sx q[0];
rz(-2.4099775) q[0];
sx q[0];
rz(-2.0425955) q[0];
rz(-pi) q[1];
x q[1];
rz(1.556384) q[2];
sx q[2];
rz(-1.3098798) q[2];
sx q[2];
rz(-0.018785611) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1085775) q[1];
sx q[1];
rz(-2.3653154) q[1];
sx q[1];
rz(-0.099221512) q[1];
rz(-0.57119675) q[3];
sx q[3];
rz(-2.9112175) q[3];
sx q[3];
rz(-0.26709589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15829076) q[2];
sx q[2];
rz(-3.0088708) q[2];
sx q[2];
rz(2.1251202) q[2];
rz(0.086183444) q[3];
sx q[3];
rz(-2.1250171) q[3];
sx q[3];
rz(0.76519722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30855274) q[0];
sx q[0];
rz(-0.85423952) q[0];
sx q[0];
rz(2.7225851) q[0];
rz(2.6047756) q[1];
sx q[1];
rz(-0.62428004) q[1];
sx q[1];
rz(0.73582617) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3296649) q[0];
sx q[0];
rz(-2.4617534) q[0];
sx q[0];
rz(0.30888866) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4707253) q[2];
sx q[2];
rz(-2.3773674) q[2];
sx q[2];
rz(-0.66167458) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.23259737) q[1];
sx q[1];
rz(-0.67282644) q[1];
sx q[1];
rz(2.7907967) q[1];
rz(-0.60634585) q[3];
sx q[3];
rz(-2.5744573) q[3];
sx q[3];
rz(-0.66886574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1824823) q[2];
sx q[2];
rz(-2.3190494) q[2];
sx q[2];
rz(0.62270069) q[2];
rz(-2.4032118) q[3];
sx q[3];
rz(-1.1850971) q[3];
sx q[3];
rz(0.27899376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5644792) q[0];
sx q[0];
rz(-0.27272419) q[0];
sx q[0];
rz(0.96584366) q[0];
rz(-2.4979757) q[1];
sx q[1];
rz(-1.8543961) q[1];
sx q[1];
rz(-2.054945) q[1];
rz(2.9109091) q[2];
sx q[2];
rz(-1.9311957) q[2];
sx q[2];
rz(0.25372505) q[2];
rz(0.52538659) q[3];
sx q[3];
rz(-2.864884) q[3];
sx q[3];
rz(-1.3340193) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
