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
rz(4.3972754) q[0];
sx q[0];
rz(9.7527405) q[0];
rz(2.9070931) q[1];
sx q[1];
rz(-0.20107888) q[1];
sx q[1];
rz(-3.0501563) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56337315) q[0];
sx q[0];
rz(-1.9679929) q[0];
sx q[0];
rz(-2.5114775) q[0];
x q[1];
rz(2.6684127) q[2];
sx q[2];
rz(-2.8521529) q[2];
sx q[2];
rz(-0.48130408) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5603197) q[1];
sx q[1];
rz(-2.8898507) q[1];
sx q[1];
rz(-1.2341577) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9760194) q[3];
sx q[3];
rz(-1.752749) q[3];
sx q[3];
rz(0.44427696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4771007) q[2];
sx q[2];
rz(-0.97350073) q[2];
sx q[2];
rz(1.1260024) q[2];
rz(-0.27515718) q[3];
sx q[3];
rz(-0.61029172) q[3];
sx q[3];
rz(-0.90308213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-3.1317516) q[0];
sx q[0];
rz(-2.4480902) q[0];
sx q[0];
rz(2.4480208) q[0];
rz(1.0961078) q[1];
sx q[1];
rz(-2.1577436) q[1];
sx q[1];
rz(-0.19031659) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0424461) q[0];
sx q[0];
rz(-2.6876246) q[0];
sx q[0];
rz(-2.0367665) q[0];
rz(2.6066577) q[2];
sx q[2];
rz(-1.2282279) q[2];
sx q[2];
rz(-1.554622) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7029876) q[1];
sx q[1];
rz(-2.4191796) q[1];
sx q[1];
rz(-1.6045251) q[1];
rz(-pi) q[2];
rz(2.4137647) q[3];
sx q[3];
rz(-1.2065294) q[3];
sx q[3];
rz(1.4621853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.50743121) q[2];
sx q[2];
rz(-2.4298411) q[2];
sx q[2];
rz(-1.921839) q[2];
rz(0.1427342) q[3];
sx q[3];
rz(-2.1295363) q[3];
sx q[3];
rz(-0.90863168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.398657) q[0];
sx q[0];
rz(-1.0899028) q[0];
sx q[0];
rz(2.3213342) q[0];
rz(-2.8495158) q[1];
sx q[1];
rz(-2.0675817) q[1];
sx q[1];
rz(1.8935727) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6662212) q[0];
sx q[0];
rz(-2.1682122) q[0];
sx q[0];
rz(-1.4587547) q[0];
x q[1];
rz(1.6349995) q[2];
sx q[2];
rz(-1.3639796) q[2];
sx q[2];
rz(-1.1484255) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.65590811) q[1];
sx q[1];
rz(-0.67042065) q[1];
sx q[1];
rz(0.1078492) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1850584) q[3];
sx q[3];
rz(-0.68813656) q[3];
sx q[3];
rz(-3.0582173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.32039207) q[2];
sx q[2];
rz(-1.0504861) q[2];
sx q[2];
rz(-1.9281663) q[2];
rz(0.16472566) q[3];
sx q[3];
rz(-2.2475593) q[3];
sx q[3];
rz(-3.059982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7320025) q[0];
sx q[0];
rz(-1.2816757) q[0];
sx q[0];
rz(2.9555292) q[0];
rz(2.9371254) q[1];
sx q[1];
rz(-0.47195131) q[1];
sx q[1];
rz(1.2971372) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9545249) q[0];
sx q[0];
rz(-1.6141119) q[0];
sx q[0];
rz(1.4403507) q[0];
rz(-1.2816216) q[2];
sx q[2];
rz(-1.4099858) q[2];
sx q[2];
rz(2.5922054) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2343826) q[1];
sx q[1];
rz(-1.0115336) q[1];
sx q[1];
rz(-2.2940966) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1277539) q[3];
sx q[3];
rz(-1.3912364) q[3];
sx q[3];
rz(1.9765215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.90594784) q[2];
sx q[2];
rz(-0.79259688) q[2];
sx q[2];
rz(-0.91919351) q[2];
rz(2.8202608) q[3];
sx q[3];
rz(-2.0733757) q[3];
sx q[3];
rz(-1.2478158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9677251) q[0];
sx q[0];
rz(-0.50650948) q[0];
sx q[0];
rz(0.87819535) q[0];
rz(1.325266) q[1];
sx q[1];
rz(-1.4122496) q[1];
sx q[1];
rz(1.426288) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8844937) q[0];
sx q[0];
rz(-1.4887267) q[0];
sx q[0];
rz(0.9135855) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8833141) q[2];
sx q[2];
rz(-1.7973571) q[2];
sx q[2];
rz(-2.8001919) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.95477415) q[1];
sx q[1];
rz(-0.84328077) q[1];
sx q[1];
rz(0.14528841) q[1];
rz(2.812643) q[3];
sx q[3];
rz(-2.3429686) q[3];
sx q[3];
rz(2.9165099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8695801) q[2];
sx q[2];
rz(-11*pi/13) q[2];
sx q[2];
rz(0.041794725) q[2];
rz(0.061491866) q[3];
sx q[3];
rz(-2.0525335) q[3];
sx q[3];
rz(0.4367691) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7866661) q[0];
sx q[0];
rz(-0.085952856) q[0];
sx q[0];
rz(2.1110995) q[0];
rz(0.73973918) q[1];
sx q[1];
rz(-1.5286427) q[1];
sx q[1];
rz(2.5700263) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.556658) q[0];
sx q[0];
rz(-2.1437763) q[0];
sx q[0];
rz(-0.95227382) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7320485) q[2];
sx q[2];
rz(-0.67220062) q[2];
sx q[2];
rz(1.6598998) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.73733444) q[1];
sx q[1];
rz(-1.3306276) q[1];
sx q[1];
rz(-2.0391697) q[1];
rz(-pi) q[2];
rz(-2.3267641) q[3];
sx q[3];
rz(-2.1216672) q[3];
sx q[3];
rz(1.8967472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.17343865) q[2];
sx q[2];
rz(-2.3996694) q[2];
sx q[2];
rz(2.5773933) q[2];
rz(0.12600222) q[3];
sx q[3];
rz(-1.6832451) q[3];
sx q[3];
rz(2.8469767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0903704) q[0];
sx q[0];
rz(-2.9416961) q[0];
sx q[0];
rz(-2.4293161) q[0];
rz(-0.5258711) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(-2.4760822) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77008477) q[0];
sx q[0];
rz(-1.5697644) q[0];
sx q[0];
rz(-1.3285711) q[0];
rz(0.87343563) q[2];
sx q[2];
rz(-1.4566112) q[2];
sx q[2];
rz(-0.21542491) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.29384089) q[1];
sx q[1];
rz(-1.7898702) q[1];
sx q[1];
rz(1.1024544) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0379167) q[3];
sx q[3];
rz(-2.8623192) q[3];
sx q[3];
rz(-2.6168407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9843288) q[2];
sx q[2];
rz(-2.4815346) q[2];
sx q[2];
rz(-1.4228014) q[2];
rz(-3.026399) q[3];
sx q[3];
rz(-2.6264103) q[3];
sx q[3];
rz(-0.19259024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0916864) q[0];
sx q[0];
rz(-2.005907) q[0];
sx q[0];
rz(-2.9507622) q[0];
rz(-2.514839) q[1];
sx q[1];
rz(-1.0160867) q[1];
sx q[1];
rz(-0.33871067) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3074293) q[0];
sx q[0];
rz(-1.2273664) q[0];
sx q[0];
rz(2.9817581) q[0];
rz(0.4523925) q[2];
sx q[2];
rz(-0.49491844) q[2];
sx q[2];
rz(2.6099043) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.45704493) q[1];
sx q[1];
rz(-2.0431879) q[1];
sx q[1];
rz(-0.12462516) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2440153) q[3];
sx q[3];
rz(-1.8339001) q[3];
sx q[3];
rz(2.9941878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4954341) q[2];
sx q[2];
rz(-0.63445264) q[2];
sx q[2];
rz(-0.70043606) q[2];
rz(-0.8979848) q[3];
sx q[3];
rz(-1.8549517) q[3];
sx q[3];
rz(-0.17351304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
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
rz(-0.53196466) q[0];
sx q[0];
rz(-0.48159972) q[0];
sx q[0];
rz(0.65834808) q[0];
rz(0.61093962) q[1];
sx q[1];
rz(-1.8860203) q[1];
sx q[1];
rz(0.13959612) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0848207) q[0];
sx q[0];
rz(-0.064282566) q[0];
sx q[0];
rz(-1.8875185) q[0];
rz(-1.6330958) q[2];
sx q[2];
rz(-0.75714105) q[2];
sx q[2];
rz(1.8119259) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3176206) q[1];
sx q[1];
rz(-2.6156153) q[1];
sx q[1];
rz(3.0082263) q[1];
rz(-pi) q[2];
rz(3.0890907) q[3];
sx q[3];
rz(-2.2146261) q[3];
sx q[3];
rz(2.1569463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5247941) q[2];
sx q[2];
rz(-2.1308265) q[2];
sx q[2];
rz(-2.810478) q[2];
rz(-0.75774276) q[3];
sx q[3];
rz(-2.7527633) q[3];
sx q[3];
rz(-0.087879114) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8452633) q[0];
sx q[0];
rz(-2.1305278) q[0];
sx q[0];
rz(-2.9560126) q[0];
rz(-2.045385) q[1];
sx q[1];
rz(-0.21462333) q[1];
sx q[1];
rz(1.4846444) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97666868) q[0];
sx q[0];
rz(-2.6161368) q[0];
sx q[0];
rz(2.1303961) q[0];
rz(-1.8469641) q[2];
sx q[2];
rz(-1.7880926) q[2];
sx q[2];
rz(1.5751788) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6757322) q[1];
sx q[1];
rz(-1.670174) q[1];
sx q[1];
rz(1.194721) q[1];
x q[2];
rz(1.5766034) q[3];
sx q[3];
rz(-2.4366637) q[3];
sx q[3];
rz(-2.4095636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.93402702) q[2];
sx q[2];
rz(-2.2521314) q[2];
sx q[2];
rz(2.5893842) q[2];
rz(-2.3637555) q[3];
sx q[3];
rz(-2.2838897) q[3];
sx q[3];
rz(-2.623693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8778397) q[0];
sx q[0];
rz(-1.629566) q[0];
sx q[0];
rz(-1.4019479) q[0];
rz(2.9539625) q[1];
sx q[1];
rz(-1.3627121) q[1];
sx q[1];
rz(-0.77308853) q[1];
rz(0.4734584) q[2];
sx q[2];
rz(-0.31242328) q[2];
sx q[2];
rz(-1.7996126) q[2];
rz(-0.68708146) q[3];
sx q[3];
rz(-1.0071181) q[3];
sx q[3];
rz(-1.3345171) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];