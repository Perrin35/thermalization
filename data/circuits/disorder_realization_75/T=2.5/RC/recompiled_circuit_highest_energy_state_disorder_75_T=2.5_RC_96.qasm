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
rz(1.6504352) q[0];
sx q[0];
rz(-0.37547922) q[0];
sx q[0];
rz(-0.39146358) q[0];
rz(0.23671167) q[1];
sx q[1];
rz(-2.1032636) q[1];
sx q[1];
rz(-2.2957323) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0856871) q[0];
sx q[0];
rz(-0.90856594) q[0];
sx q[0];
rz(-1.7895749) q[0];
rz(0.45028668) q[2];
sx q[2];
rz(-1.1495207) q[2];
sx q[2];
rz(-2.7305528) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9917984) q[1];
sx q[1];
rz(-1.845896) q[1];
sx q[1];
rz(-2.8582948) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3491783) q[3];
sx q[3];
rz(-1.0842348) q[3];
sx q[3];
rz(0.38565991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.55964959) q[2];
sx q[2];
rz(-1.0407642) q[2];
sx q[2];
rz(-2.9259658) q[2];
rz(-0.98254472) q[3];
sx q[3];
rz(-1.6590786) q[3];
sx q[3];
rz(-1.2561579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7243645) q[0];
sx q[0];
rz(-0.094014458) q[0];
sx q[0];
rz(-2.7650058) q[0];
rz(-1.6603893) q[1];
sx q[1];
rz(-1.9638289) q[1];
sx q[1];
rz(-2.1362163) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.004285) q[0];
sx q[0];
rz(-1.2457324) q[0];
sx q[0];
rz(0.0079197366) q[0];
rz(-pi) q[1];
rz(-0.89400228) q[2];
sx q[2];
rz(-1.4419705) q[2];
sx q[2];
rz(2.1137546) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.59170106) q[1];
sx q[1];
rz(-1.5366067) q[1];
sx q[1];
rz(-0.66877613) q[1];
rz(-pi) q[2];
rz(-1.1930258) q[3];
sx q[3];
rz(-2.0888512) q[3];
sx q[3];
rz(-2.4729721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5542095) q[2];
sx q[2];
rz(-2.6770834) q[2];
sx q[2];
rz(-1.0760388) q[2];
rz(0.46766034) q[3];
sx q[3];
rz(-0.83100072) q[3];
sx q[3];
rz(0.20461288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49863807) q[0];
sx q[0];
rz(-1.0774281) q[0];
sx q[0];
rz(2.0306008) q[0];
rz(0.30771646) q[1];
sx q[1];
rz(-1.6650763) q[1];
sx q[1];
rz(1.3002546) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0205524) q[0];
sx q[0];
rz(-2.5835134) q[0];
sx q[0];
rz(0.16772018) q[0];
x q[1];
rz(2.563649) q[2];
sx q[2];
rz(-1.5405077) q[2];
sx q[2];
rz(2.2453215) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.96836263) q[1];
sx q[1];
rz(-0.97665962) q[1];
sx q[1];
rz(0.51715516) q[1];
x q[2];
rz(-2.6407317) q[3];
sx q[3];
rz(-1.259049) q[3];
sx q[3];
rz(-0.65697981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51851455) q[2];
sx q[2];
rz(-2.0936421) q[2];
sx q[2];
rz(0.19185647) q[2];
rz(-2.0375371) q[3];
sx q[3];
rz(-2.7147229) q[3];
sx q[3];
rz(-1.2087315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0693102) q[0];
sx q[0];
rz(-2.5492714) q[0];
sx q[0];
rz(2.2533672) q[0];
rz(-2.8690673) q[1];
sx q[1];
rz(-1.2963632) q[1];
sx q[1];
rz(-1.1558418) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4711535) q[0];
sx q[0];
rz(-1.5799755) q[0];
sx q[0];
rz(1.8646452) q[0];
x q[1];
rz(0.49057284) q[2];
sx q[2];
rz(-1.6075879) q[2];
sx q[2];
rz(1.8951777) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.19570505) q[1];
sx q[1];
rz(-0.81521704) q[1];
sx q[1];
rz(2.7570711) q[1];
rz(2.6937595) q[3];
sx q[3];
rz(-1.5671697) q[3];
sx q[3];
rz(-1.8254395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.83879519) q[2];
sx q[2];
rz(-0.28442997) q[2];
sx q[2];
rz(-2.3885041) q[2];
rz(-2.6436842) q[3];
sx q[3];
rz(-1.6585191) q[3];
sx q[3];
rz(-0.30043093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4860151) q[0];
sx q[0];
rz(-0.97633728) q[0];
sx q[0];
rz(-1.8560386) q[0];
rz(-1.1898419) q[1];
sx q[1];
rz(-2.4557476) q[1];
sx q[1];
rz(1.5600342) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6046977) q[0];
sx q[0];
rz(-1.6486537) q[0];
sx q[0];
rz(-1.4417062) q[0];
rz(2.6597775) q[2];
sx q[2];
rz(-0.31793943) q[2];
sx q[2];
rz(-0.16127333) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6717357) q[1];
sx q[1];
rz(-1.7845881) q[1];
sx q[1];
rz(-0.44624568) q[1];
x q[2];
rz(-0.26435477) q[3];
sx q[3];
rz(-1.3761545) q[3];
sx q[3];
rz(-2.1771012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.84216422) q[2];
sx q[2];
rz(-0.93141586) q[2];
sx q[2];
rz(0.13775873) q[2];
rz(-0.44451851) q[3];
sx q[3];
rz(-1.3714182) q[3];
sx q[3];
rz(-0.81454149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-2.923362) q[0];
sx q[0];
rz(-2.2683999) q[0];
sx q[0];
rz(2.5905304) q[0];
rz(2.8463544) q[1];
sx q[1];
rz(-2.4199838) q[1];
sx q[1];
rz(-0.74337426) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6310204) q[0];
sx q[0];
rz(-1.6845682) q[0];
sx q[0];
rz(-1.7049417) q[0];
x q[1];
rz(-0.34140519) q[2];
sx q[2];
rz(-1.3479173) q[2];
sx q[2];
rz(2.6472732) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.049014) q[1];
sx q[1];
rz(-0.72737073) q[1];
sx q[1];
rz(-1.7197827) q[1];
x q[2];
rz(-2.4427209) q[3];
sx q[3];
rz(-1.0632319) q[3];
sx q[3];
rz(-2.5273049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.729852) q[2];
sx q[2];
rz(-0.19062947) q[2];
sx q[2];
rz(2.8594678) q[2];
rz(2.7052346) q[3];
sx q[3];
rz(-1.8264419) q[3];
sx q[3];
rz(0.25562975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.2082763) q[0];
sx q[0];
rz(-2.0774807) q[0];
sx q[0];
rz(2.5489885) q[0];
rz(-2.0049877) q[1];
sx q[1];
rz(-1.2717783) q[1];
sx q[1];
rz(1.0677387) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16780218) q[0];
sx q[0];
rz(-2.4616957) q[0];
sx q[0];
rz(2.8907084) q[0];
rz(-2.7726658) q[2];
sx q[2];
rz(-1.6875982) q[2];
sx q[2];
rz(2.359085) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.021879747) q[1];
sx q[1];
rz(-1.2924866) q[1];
sx q[1];
rz(1.0779508) q[1];
rz(-pi) q[2];
rz(2.7619902) q[3];
sx q[3];
rz(-2.6983454) q[3];
sx q[3];
rz(-0.28397181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9390255) q[2];
sx q[2];
rz(-1.2747526) q[2];
sx q[2];
rz(2.1620046) q[2];
rz(-3.0136287) q[3];
sx q[3];
rz(-2.1991859) q[3];
sx q[3];
rz(-1.6894107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1487883) q[0];
sx q[0];
rz(-0.95071852) q[0];
sx q[0];
rz(0.08091452) q[0];
rz(3.0203536) q[1];
sx q[1];
rz(-2.464005) q[1];
sx q[1];
rz(-2.0176719) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1428536) q[0];
sx q[0];
rz(-0.8585728) q[0];
sx q[0];
rz(-1.0741878) q[0];
rz(-pi) q[1];
rz(-0.10350714) q[2];
sx q[2];
rz(-2.7324504) q[2];
sx q[2];
rz(-1.658866) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7940136) q[1];
sx q[1];
rz(-1.4718143) q[1];
sx q[1];
rz(-0.35225676) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2865661) q[3];
sx q[3];
rz(-2.5919302) q[3];
sx q[3];
rz(-0.46297234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8671882) q[2];
sx q[2];
rz(-0.97325456) q[2];
sx q[2];
rz(1.6483866) q[2];
rz(-0.46227208) q[3];
sx q[3];
rz(-0.8148163) q[3];
sx q[3];
rz(-1.7260684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3129775) q[0];
sx q[0];
rz(-1.7531489) q[0];
sx q[0];
rz(1.9762565) q[0];
rz(-0.041821592) q[1];
sx q[1];
rz(-2.1526497) q[1];
sx q[1];
rz(-1.8080541) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8671252) q[0];
sx q[0];
rz(-2.2750912) q[0];
sx q[0];
rz(-0.88732669) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5749349) q[2];
sx q[2];
rz(-1.88961) q[2];
sx q[2];
rz(-0.37296527) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.11631498) q[1];
sx q[1];
rz(-1.9851369) q[1];
sx q[1];
rz(-1.6700909) q[1];
rz(0.72527146) q[3];
sx q[3];
rz(-0.15003157) q[3];
sx q[3];
rz(-1.6440441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8970783) q[2];
sx q[2];
rz(-1.8316734) q[2];
sx q[2];
rz(-0.46373996) q[2];
rz(2.5044299) q[3];
sx q[3];
rz(-1.623268) q[3];
sx q[3];
rz(0.37850982) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8265182) q[0];
sx q[0];
rz(-1.0239064) q[0];
sx q[0];
rz(-1.4950604) q[0];
rz(0.82978326) q[1];
sx q[1];
rz(-1.8779495) q[1];
sx q[1];
rz(1.319938) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0129725) q[0];
sx q[0];
rz(-1.9397598) q[0];
sx q[0];
rz(1.7075383) q[0];
rz(-pi) q[1];
rz(2.7097335) q[2];
sx q[2];
rz(-1.7786025) q[2];
sx q[2];
rz(2.6450405) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.9480919) q[1];
sx q[1];
rz(-1.337916) q[1];
sx q[1];
rz(-1.3629254) q[1];
rz(1.9687555) q[3];
sx q[3];
rz(-0.59919849) q[3];
sx q[3];
rz(0.20922449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39184555) q[2];
sx q[2];
rz(-1.7228935) q[2];
sx q[2];
rz(-2.6302272) q[2];
rz(0.5558719) q[3];
sx q[3];
rz(-1.116773) q[3];
sx q[3];
rz(2.6218124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2340672) q[0];
sx q[0];
rz(-3.0341442) q[0];
sx q[0];
rz(-2.7035614) q[0];
rz(3.0829433) q[1];
sx q[1];
rz(-1.6390683) q[1];
sx q[1];
rz(-1.0674089) q[1];
rz(3.0728523) q[2];
sx q[2];
rz(-0.70211997) q[2];
sx q[2];
rz(1.4236249) q[2];
rz(-2.7300446) q[3];
sx q[3];
rz(-0.66933142) q[3];
sx q[3];
rz(0.29585024) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
