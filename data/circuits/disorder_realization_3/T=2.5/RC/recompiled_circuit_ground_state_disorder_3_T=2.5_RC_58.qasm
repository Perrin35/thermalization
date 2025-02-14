OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8752911) q[0];
sx q[0];
rz(-2.7346791) q[0];
sx q[0];
rz(-0.20242515) q[0];
rz(1.8316733) q[1];
sx q[1];
rz(4.2869422) q[1];
sx q[1];
rz(11.577679) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3920445) q[0];
sx q[0];
rz(-0.77717669) q[0];
sx q[0];
rz(-1.8424319) q[0];
rz(-pi) q[1];
rz(1.7086171) q[2];
sx q[2];
rz(-1.4933961) q[2];
sx q[2];
rz(-3.1397506) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1867849) q[1];
sx q[1];
rz(-2.4802172) q[1];
sx q[1];
rz(-1.1977045) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0172431) q[3];
sx q[3];
rz(-1.9938339) q[3];
sx q[3];
rz(1.7250586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.53628355) q[2];
sx q[2];
rz(-2.3399957) q[2];
sx q[2];
rz(2.0901704) q[2];
rz(-1.5031523) q[3];
sx q[3];
rz(-1.7283311) q[3];
sx q[3];
rz(0.2352636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94238344) q[0];
sx q[0];
rz(-1.0301882) q[0];
sx q[0];
rz(0.27401608) q[0];
rz(2.7577397) q[1];
sx q[1];
rz(-0.63261384) q[1];
sx q[1];
rz(-1.8208108) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1766168) q[0];
sx q[0];
rz(-2.9701485) q[0];
sx q[0];
rz(1.9587396) q[0];
rz(-pi) q[1];
rz(2.6422149) q[2];
sx q[2];
rz(-1.3670397) q[2];
sx q[2];
rz(2.8600542) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.40055028) q[1];
sx q[1];
rz(-2.2930751) q[1];
sx q[1];
rz(-3.0781915) q[1];
rz(-pi) q[2];
rz(-2.3736773) q[3];
sx q[3];
rz(-0.95988279) q[3];
sx q[3];
rz(-3.0141941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.5968093) q[2];
sx q[2];
rz(-1.7388672) q[2];
sx q[2];
rz(1.7421494) q[2];
rz(2.5341189) q[3];
sx q[3];
rz(-1.4676658) q[3];
sx q[3];
rz(-2.7510711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6121599) q[0];
sx q[0];
rz(-2.2859892) q[0];
sx q[0];
rz(2.8223619) q[0];
rz(0.051479738) q[1];
sx q[1];
rz(-1.9620644) q[1];
sx q[1];
rz(2.0850339) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7565472) q[0];
sx q[0];
rz(-2.9699885) q[0];
sx q[0];
rz(1.8770939) q[0];
x q[1];
rz(-2.276307) q[2];
sx q[2];
rz(-1.6076536) q[2];
sx q[2];
rz(-1.8344022) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.776538) q[1];
sx q[1];
rz(-2.1927196) q[1];
sx q[1];
rz(-2.9261241) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7848739) q[3];
sx q[3];
rz(-0.41164648) q[3];
sx q[3];
rz(-2.1901707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1928548) q[2];
sx q[2];
rz(-1.4317908) q[2];
sx q[2];
rz(-1.3445492) q[2];
rz(-2.2236842) q[3];
sx q[3];
rz(-2.5036006) q[3];
sx q[3];
rz(1.8857672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6767839) q[0];
sx q[0];
rz(-2.5510241) q[0];
sx q[0];
rz(1.6326686) q[0];
rz(-2.6347939) q[1];
sx q[1];
rz(-1.9281887) q[1];
sx q[1];
rz(-0.39055821) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4669663) q[0];
sx q[0];
rz(-1.5780562) q[0];
sx q[0];
rz(1.566559) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3803595) q[2];
sx q[2];
rz(-1.197063) q[2];
sx q[2];
rz(-0.98183364) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4949172) q[1];
sx q[1];
rz(-0.61958379) q[1];
sx q[1];
rz(-1.619672) q[1];
rz(-pi) q[2];
rz(-0.5162913) q[3];
sx q[3];
rz(-2.9598979) q[3];
sx q[3];
rz(-2.5207375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.249923) q[2];
sx q[2];
rz(-0.93081433) q[2];
sx q[2];
rz(2.2451952) q[2];
rz(-2.2253288) q[3];
sx q[3];
rz(-0.93583411) q[3];
sx q[3];
rz(0.96074444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59178281) q[0];
sx q[0];
rz(-0.35583219) q[0];
sx q[0];
rz(0.47019666) q[0];
rz(1.7377986) q[1];
sx q[1];
rz(-2.4411968) q[1];
sx q[1];
rz(1.2276924) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5915697) q[0];
sx q[0];
rz(-1.8409022) q[0];
sx q[0];
rz(3.1349934) q[0];
rz(-pi) q[1];
rz(0.63666261) q[2];
sx q[2];
rz(-2.4877254) q[2];
sx q[2];
rz(0.9034397) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2254082) q[1];
sx q[1];
rz(-2.492444) q[1];
sx q[1];
rz(-1.445392) q[1];
rz(-pi) q[2];
rz(-1.8813666) q[3];
sx q[3];
rz(-1.0538007) q[3];
sx q[3];
rz(2.9152512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.97051835) q[2];
sx q[2];
rz(-0.71275622) q[2];
sx q[2];
rz(1.3842545) q[2];
rz(-0.61128831) q[3];
sx q[3];
rz(-1.7620112) q[3];
sx q[3];
rz(-2.7454564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1432081) q[0];
sx q[0];
rz(-1.0118326) q[0];
sx q[0];
rz(0.55214733) q[0];
rz(0.72969189) q[1];
sx q[1];
rz(-2.153986) q[1];
sx q[1];
rz(-2.139835) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14267966) q[0];
sx q[0];
rz(-0.15593869) q[0];
sx q[0];
rz(-2.1151311) q[0];
rz(-pi) q[1];
rz(-0.4977141) q[2];
sx q[2];
rz(-0.98955284) q[2];
sx q[2];
rz(-2.0037985) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3730091) q[1];
sx q[1];
rz(-0.78861744) q[1];
sx q[1];
rz(2.8810701) q[1];
rz(-pi) q[2];
rz(-2.5252063) q[3];
sx q[3];
rz(-2.3616543) q[3];
sx q[3];
rz(-2.9857948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.45727649) q[2];
sx q[2];
rz(-1.7924954) q[2];
sx q[2];
rz(-1.7982193) q[2];
rz(0.7243048) q[3];
sx q[3];
rz(-1.2035921) q[3];
sx q[3];
rz(0.30266416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-2.2398719) q[0];
sx q[0];
rz(-1.3322823) q[0];
sx q[0];
rz(0.72541642) q[0];
rz(1.8431009) q[1];
sx q[1];
rz(-0.45227554) q[1];
sx q[1];
rz(2.8628912) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50510079) q[0];
sx q[0];
rz(-2.1064225) q[0];
sx q[0];
rz(-1.7477751) q[0];
x q[1];
rz(-2.8093074) q[2];
sx q[2];
rz(-1.1934308) q[2];
sx q[2];
rz(-1.8695631) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.78165752) q[1];
sx q[1];
rz(-1.2664478) q[1];
sx q[1];
rz(-2.3001647) q[1];
rz(-pi) q[2];
rz(-0.82741957) q[3];
sx q[3];
rz(-1.9636834) q[3];
sx q[3];
rz(-3.06638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.38410386) q[2];
sx q[2];
rz(-0.266738) q[2];
sx q[2];
rz(0.98814803) q[2];
rz(-0.10495505) q[3];
sx q[3];
rz(-1.2982488) q[3];
sx q[3];
rz(2.9102563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4526378) q[0];
sx q[0];
rz(-0.80113688) q[0];
sx q[0];
rz(2.5222006) q[0];
rz(-3.1267005) q[1];
sx q[1];
rz(-2.2512524) q[1];
sx q[1];
rz(-2.4598918) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0571201) q[0];
sx q[0];
rz(-2.1051877) q[0];
sx q[0];
rz(0.33584331) q[0];
rz(-1.6729082) q[2];
sx q[2];
rz(-2.4739304) q[2];
sx q[2];
rz(0.84102902) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.926103) q[1];
sx q[1];
rz(-1.5640096) q[1];
sx q[1];
rz(2.7899963) q[1];
x q[2];
rz(2.9497164) q[3];
sx q[3];
rz(-0.8423281) q[3];
sx q[3];
rz(1.6863562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2316078) q[2];
sx q[2];
rz(-2.3384428) q[2];
sx q[2];
rz(-0.33965674) q[2];
rz(-0.61399442) q[3];
sx q[3];
rz(-1.3795779) q[3];
sx q[3];
rz(1.5617255) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9851819) q[0];
sx q[0];
rz(-2.7934533) q[0];
sx q[0];
rz(2.2871995) q[0];
rz(-0.59335452) q[1];
sx q[1];
rz(-2.1095468) q[1];
sx q[1];
rz(2.8462483) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90705425) q[0];
sx q[0];
rz(-3.0772644) q[0];
sx q[0];
rz(2.6049019) q[0];
rz(-pi) q[1];
rz(-1.5200843) q[2];
sx q[2];
rz(-0.4182294) q[2];
sx q[2];
rz(2.744368) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.015739) q[1];
sx q[1];
rz(-2.2496114) q[1];
sx q[1];
rz(-2.5126481) q[1];
rz(-1.2611102) q[3];
sx q[3];
rz(-1.4638374) q[3];
sx q[3];
rz(-0.84236077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.16193834) q[2];
sx q[2];
rz(-1.4583959) q[2];
sx q[2];
rz(2.4597607) q[2];
rz(1.2006867) q[3];
sx q[3];
rz(-2.3537945) q[3];
sx q[3];
rz(0.43584263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1107263) q[0];
sx q[0];
rz(-1.2799355) q[0];
sx q[0];
rz(-0.50586343) q[0];
rz(1.0181631) q[1];
sx q[1];
rz(-1.70645) q[1];
sx q[1];
rz(0.86046576) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9415362) q[0];
sx q[0];
rz(-2.464009) q[0];
sx q[0];
rz(-0.075417544) q[0];
rz(-pi) q[1];
rz(-1.127709) q[2];
sx q[2];
rz(-2.8193316) q[2];
sx q[2];
rz(1.6784422) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.082749359) q[1];
sx q[1];
rz(-0.25139055) q[1];
sx q[1];
rz(-0.13240929) q[1];
rz(2.5768336) q[3];
sx q[3];
rz(-1.0339206) q[3];
sx q[3];
rz(2.2220461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0818417) q[2];
sx q[2];
rz(-1.1885208) q[2];
sx q[2];
rz(-0.58319485) q[2];
rz(3.1405295) q[3];
sx q[3];
rz(-0.93102264) q[3];
sx q[3];
rz(-2.6482436) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6625593) q[0];
sx q[0];
rz(-1.722535) q[0];
sx q[0];
rz(-2.2395635) q[0];
rz(2.5262911) q[1];
sx q[1];
rz(-1.4882144) q[1];
sx q[1];
rz(-1.8019713) q[1];
rz(-0.42846305) q[2];
sx q[2];
rz(-0.56461538) q[2];
sx q[2];
rz(2.3587894) q[2];
rz(-2.9499182) q[3];
sx q[3];
rz(-2.1379466) q[3];
sx q[3];
rz(2.9320075) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
