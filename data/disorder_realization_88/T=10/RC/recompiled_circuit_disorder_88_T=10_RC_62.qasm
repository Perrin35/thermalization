OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.24703439) q[0];
sx q[0];
rz(-1.8859099) q[0];
sx q[0];
rz(0.32796252) q[0];
rz(-0.2344996) q[1];
sx q[1];
rz(3.3426715) q[1];
sx q[1];
rz(9.3333416) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56337315) q[0];
sx q[0];
rz(-1.9679929) q[0];
sx q[0];
rz(2.5114775) q[0];
x q[1];
rz(-2.6684127) q[2];
sx q[2];
rz(-0.28943974) q[2];
sx q[2];
rz(2.6602886) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5603197) q[1];
sx q[1];
rz(-2.8898507) q[1];
sx q[1];
rz(1.2341577) q[1];
rz(-pi) q[2];
x q[2];
rz(2.923008) q[3];
sx q[3];
rz(-0.9871452) q[3];
sx q[3];
rz(-2.1368795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.66449195) q[2];
sx q[2];
rz(-2.1680919) q[2];
sx q[2];
rz(-2.0155902) q[2];
rz(0.27515718) q[3];
sx q[3];
rz(-0.61029172) q[3];
sx q[3];
rz(0.90308213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1317516) q[0];
sx q[0];
rz(-0.69350243) q[0];
sx q[0];
rz(2.4480208) q[0];
rz(-2.0454848) q[1];
sx q[1];
rz(-0.98384905) q[1];
sx q[1];
rz(-2.9512761) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0424461) q[0];
sx q[0];
rz(-2.6876246) q[0];
sx q[0];
rz(2.0367665) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1778189) q[2];
sx q[2];
rz(-2.071638) q[2];
sx q[2];
rz(-2.9289392) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.43860501) q[1];
sx q[1];
rz(-2.4191796) q[1];
sx q[1];
rz(-1.6045251) q[1];
rz(-0.72782794) q[3];
sx q[3];
rz(-1.2065294) q[3];
sx q[3];
rz(-1.6794074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6341614) q[2];
sx q[2];
rz(-0.71175152) q[2];
sx q[2];
rz(-1.2197536) q[2];
rz(-2.9988585) q[3];
sx q[3];
rz(-2.1295363) q[3];
sx q[3];
rz(2.232961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-2.398657) q[0];
sx q[0];
rz(-1.0899028) q[0];
sx q[0];
rz(2.3213342) q[0];
rz(2.8495158) q[1];
sx q[1];
rz(-1.074011) q[1];
sx q[1];
rz(-1.2480199) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67278636) q[0];
sx q[0];
rz(-2.5350223) q[0];
sx q[0];
rz(2.9787105) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6349995) q[2];
sx q[2];
rz(-1.777613) q[2];
sx q[2];
rz(-1.1484255) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4856845) q[1];
sx q[1];
rz(-0.67042065) q[1];
sx q[1];
rz(0.1078492) q[1];
x q[2];
rz(2.8415801) q[3];
sx q[3];
rz(-0.94167751) q[3];
sx q[3];
rz(-0.56738561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8212006) q[2];
sx q[2];
rz(-2.0911066) q[2];
sx q[2];
rz(-1.9281663) q[2];
rz(2.976867) q[3];
sx q[3];
rz(-2.2475593) q[3];
sx q[3];
rz(3.059982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7320025) q[0];
sx q[0];
rz(-1.2816757) q[0];
sx q[0];
rz(2.9555292) q[0];
rz(-0.20446725) q[1];
sx q[1];
rz(-0.47195131) q[1];
sx q[1];
rz(-1.8444555) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38940934) q[0];
sx q[0];
rz(-1.4404738) q[0];
sx q[0];
rz(-0.043686314) q[0];
rz(-1.2816216) q[2];
sx q[2];
rz(-1.7316069) q[2];
sx q[2];
rz(-2.5922054) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.101673) q[1];
sx q[1];
rz(-2.1664201) q[1];
sx q[1];
rz(-2.4458829) q[1];
rz(1.6468871) q[3];
sx q[3];
rz(-2.9615059) q[3];
sx q[3];
rz(1.0877346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2356448) q[2];
sx q[2];
rz(-0.79259688) q[2];
sx q[2];
rz(-2.2223991) q[2];
rz(-2.8202608) q[3];
sx q[3];
rz(-1.0682169) q[3];
sx q[3];
rz(-1.2478158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17386757) q[0];
sx q[0];
rz(-2.6350832) q[0];
sx q[0];
rz(0.87819535) q[0];
rz(1.8163266) q[1];
sx q[1];
rz(-1.4122496) q[1];
sx q[1];
rz(-1.426288) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.250524) q[0];
sx q[0];
rz(-0.91618012) q[0];
sx q[0];
rz(3.0380681) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9039188) q[2];
sx q[2];
rz(-1.2665247) q[2];
sx q[2];
rz(-1.8397457) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6225699) q[1];
sx q[1];
rz(-1.4624603) q[1];
sx q[1];
rz(-2.3035754) q[1];
rz(-0.77107314) q[3];
sx q[3];
rz(-1.3372476) q[3];
sx q[3];
rz(-1.1119103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8695801) q[2];
sx q[2];
rz(-11*pi/13) q[2];
sx q[2];
rz(-3.0997979) q[2];
rz(-0.061491866) q[3];
sx q[3];
rz(-1.0890591) q[3];
sx q[3];
rz(0.4367691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.3549266) q[0];
sx q[0];
rz(-3.0556398) q[0];
sx q[0];
rz(-1.0304931) q[0];
rz(-0.73973918) q[1];
sx q[1];
rz(-1.5286427) q[1];
sx q[1];
rz(-2.5700263) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.556658) q[0];
sx q[0];
rz(-0.99781636) q[0];
sx q[0];
rz(2.1893188) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4095441) q[2];
sx q[2];
rz(-2.469392) q[2];
sx q[2];
rz(-1.4816928) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7476269) q[1];
sx q[1];
rz(-0.52226258) q[1];
sx q[1];
rz(2.0678492) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70116331) q[3];
sx q[3];
rz(-0.94651604) q[3];
sx q[3];
rz(0.13247709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17343865) q[2];
sx q[2];
rz(-0.74192321) q[2];
sx q[2];
rz(0.56419939) q[2];
rz(-0.12600222) q[3];
sx q[3];
rz(-1.6832451) q[3];
sx q[3];
rz(0.29461598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.0512222) q[0];
sx q[0];
rz(-0.19989656) q[0];
sx q[0];
rz(-2.4293161) q[0];
rz(2.6157216) q[1];
sx q[1];
rz(-0.41627517) q[1];
sx q[1];
rz(2.4760822) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3715079) q[0];
sx q[0];
rz(-1.5718282) q[0];
sx q[0];
rz(-1.3285711) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87343563) q[2];
sx q[2];
rz(-1.4566112) q[2];
sx q[2];
rz(-2.9261677) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1674541) q[1];
sx q[1];
rz(-2.0270837) q[1];
sx q[1];
rz(2.8970701) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1036759) q[3];
sx q[3];
rz(-2.8623192) q[3];
sx q[3];
rz(0.52475196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9843288) q[2];
sx q[2];
rz(-0.66005808) q[2];
sx q[2];
rz(-1.4228014) q[2];
rz(-3.026399) q[3];
sx q[3];
rz(-0.51518232) q[3];
sx q[3];
rz(0.19259024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0916864) q[0];
sx q[0];
rz(-2.005907) q[0];
sx q[0];
rz(2.9507622) q[0];
rz(0.62675369) q[1];
sx q[1];
rz(-1.0160867) q[1];
sx q[1];
rz(-0.33871067) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9324547) q[0];
sx q[0];
rz(-1.4203686) q[0];
sx q[0];
rz(-1.9183137) q[0];
x q[1];
rz(-1.3391101) q[2];
sx q[2];
rz(-1.1294239) q[2];
sx q[2];
rz(-2.1053134) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.45704493) q[1];
sx q[1];
rz(-1.0984048) q[1];
sx q[1];
rz(-3.0169675) q[1];
rz(1.2440153) q[3];
sx q[3];
rz(-1.8339001) q[3];
sx q[3];
rz(-0.14740482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.64615858) q[2];
sx q[2];
rz(-2.50714) q[2];
sx q[2];
rz(-2.4411566) q[2];
rz(2.2436079) q[3];
sx q[3];
rz(-1.8549517) q[3];
sx q[3];
rz(-0.17351304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53196466) q[0];
sx q[0];
rz(-2.6599929) q[0];
sx q[0];
rz(-2.4832446) q[0];
rz(-0.61093962) q[1];
sx q[1];
rz(-1.8860203) q[1];
sx q[1];
rz(-0.13959612) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0848207) q[0];
sx q[0];
rz(-3.0773101) q[0];
sx q[0];
rz(-1.8875185) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.058768674) q[2];
sx q[2];
rz(-2.3261056) q[2];
sx q[2];
rz(1.8975443) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3176206) q[1];
sx q[1];
rz(-2.6156153) q[1];
sx q[1];
rz(-0.13336639) q[1];
rz(-3.0890907) q[3];
sx q[3];
rz(-0.92696654) q[3];
sx q[3];
rz(-0.98464636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5247941) q[2];
sx q[2];
rz(-1.0107661) q[2];
sx q[2];
rz(2.810478) q[2];
rz(-0.75774276) q[3];
sx q[3];
rz(-0.38882935) q[3];
sx q[3];
rz(0.087879114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2963294) q[0];
sx q[0];
rz(-1.0110649) q[0];
sx q[0];
rz(0.18558003) q[0];
rz(-1.0962076) q[1];
sx q[1];
rz(-2.9269693) q[1];
sx q[1];
rz(1.4846444) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.164924) q[0];
sx q[0];
rz(-0.52545588) q[0];
sx q[0];
rz(-2.1303961) q[0];
rz(0.22557232) q[2];
sx q[2];
rz(-1.301287) q[2];
sx q[2];
rz(-3.0849506) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1440925) q[1];
sx q[1];
rz(-1.1966685) q[1];
sx q[1];
rz(-3.0348026) q[1];
rz(-pi) q[2];
rz(0.86587571) q[3];
sx q[3];
rz(-1.5745592) q[3];
sx q[3];
rz(-2.3072484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.93402702) q[2];
sx q[2];
rz(-2.2521314) q[2];
sx q[2];
rz(0.55220848) q[2];
rz(-0.77783716) q[3];
sx q[3];
rz(-0.85770291) q[3];
sx q[3];
rz(0.51789969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26375297) q[0];
sx q[0];
rz(-1.629566) q[0];
sx q[0];
rz(-1.4019479) q[0];
rz(2.9539625) q[1];
sx q[1];
rz(-1.3627121) q[1];
sx q[1];
rz(-0.77308853) q[1];
rz(-2.6681343) q[2];
sx q[2];
rz(-0.31242328) q[2];
sx q[2];
rz(-1.7996126) q[2];
rz(-2.3579303) q[3];
sx q[3];
rz(-2.2829934) q[3];
sx q[3];
rz(-0.34096277) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];