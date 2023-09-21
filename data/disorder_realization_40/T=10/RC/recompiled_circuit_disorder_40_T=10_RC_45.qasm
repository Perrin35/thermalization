OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6089132) q[0];
sx q[0];
rz(-0.37663868) q[0];
sx q[0];
rz(0.11178804) q[0];
rz(1.6821661) q[1];
sx q[1];
rz(4.7987727) q[1];
sx q[1];
rz(6.12943) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.789334) q[0];
sx q[0];
rz(-2.6129122) q[0];
sx q[0];
rz(0.60469158) q[0];
x q[1];
rz(-1.8741688) q[2];
sx q[2];
rz(-0.25250013) q[2];
sx q[2];
rz(-1.7054103) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44126025) q[1];
sx q[1];
rz(-1.2520737) q[1];
sx q[1];
rz(3.1096427) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8044756) q[3];
sx q[3];
rz(-1.9883336) q[3];
sx q[3];
rz(-1.9330213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6979606) q[2];
sx q[2];
rz(-1.7093095) q[2];
sx q[2];
rz(1.704818) q[2];
rz(2.4076961) q[3];
sx q[3];
rz(-1.5489483) q[3];
sx q[3];
rz(2.6255887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2519418) q[0];
sx q[0];
rz(-1.2263068) q[0];
sx q[0];
rz(0.92457986) q[0];
rz(2.1444767) q[1];
sx q[1];
rz(-2.6328502) q[1];
sx q[1];
rz(-1.8181713) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.130587) q[0];
sx q[0];
rz(-2.4936094) q[0];
sx q[0];
rz(2.381071) q[0];
rz(-pi) q[1];
rz(1.5263444) q[2];
sx q[2];
rz(-1.9200846) q[2];
sx q[2];
rz(-0.312422) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.70222774) q[1];
sx q[1];
rz(-0.80184466) q[1];
sx q[1];
rz(0.37978362) q[1];
x q[2];
rz(2.1453342) q[3];
sx q[3];
rz(-1.901445) q[3];
sx q[3];
rz(-1.4853256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9374342) q[2];
sx q[2];
rz(-1.5896475) q[2];
sx q[2];
rz(0.75817529) q[2];
rz(2.5126863) q[3];
sx q[3];
rz(-0.40142504) q[3];
sx q[3];
rz(1.1531856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.4085098) q[0];
sx q[0];
rz(-1.874431) q[0];
sx q[0];
rz(-2.4531903) q[0];
rz(-3.0738661) q[1];
sx q[1];
rz(-1.7522782) q[1];
sx q[1];
rz(0.53007954) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1353697) q[0];
sx q[0];
rz(-2.831499) q[0];
sx q[0];
rz(-1.9276516) q[0];
rz(-pi) q[1];
rz(-1.0447787) q[2];
sx q[2];
rz(-1.9472497) q[2];
sx q[2];
rz(-1.2029755) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0878151) q[1];
sx q[1];
rz(-0.74790819) q[1];
sx q[1];
rz(1.0653711) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.097966627) q[3];
sx q[3];
rz(-1.8091822) q[3];
sx q[3];
rz(-2.1408248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.80785859) q[2];
sx q[2];
rz(-3.1299751) q[2];
sx q[2];
rz(-2.2401436) q[2];
rz(-0.83550134) q[3];
sx q[3];
rz(-1.52799) q[3];
sx q[3];
rz(1.3114595) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9505342) q[0];
sx q[0];
rz(-1.598851) q[0];
sx q[0];
rz(0.78432551) q[0];
rz(-3.0803608) q[1];
sx q[1];
rz(-0.71413723) q[1];
sx q[1];
rz(0.13664666) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.559583) q[0];
sx q[0];
rz(-2.0499381) q[0];
sx q[0];
rz(0.14736202) q[0];
rz(-pi) q[1];
rz(2.5606974) q[2];
sx q[2];
rz(-0.53933203) q[2];
sx q[2];
rz(-2.5748594) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.20679684) q[1];
sx q[1];
rz(-3.0691642) q[1];
sx q[1];
rz(-1.9070684) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1468048) q[3];
sx q[3];
rz(-1.6487062) q[3];
sx q[3];
rz(-0.60914492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0121997) q[2];
sx q[2];
rz(-0.94038525) q[2];
sx q[2];
rz(2.5811035) q[2];
rz(-3.1292606) q[3];
sx q[3];
rz(-2.2380232) q[3];
sx q[3];
rz(-2.0509317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7085768) q[0];
sx q[0];
rz(-2.5362159) q[0];
sx q[0];
rz(0.82114712) q[0];
rz(0.87617809) q[1];
sx q[1];
rz(-2.2416302) q[1];
sx q[1];
rz(1.4076153) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8713184) q[0];
sx q[0];
rz(-0.5373913) q[0];
sx q[0];
rz(-2.4094765) q[0];
rz(2.8367923) q[2];
sx q[2];
rz(-1.0101724) q[2];
sx q[2];
rz(1.0501109) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.69406063) q[1];
sx q[1];
rz(-2.5563055) q[1];
sx q[1];
rz(-1.3259757) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2110662) q[3];
sx q[3];
rz(-1.3621646) q[3];
sx q[3];
rz(-1.2071351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6901107) q[2];
sx q[2];
rz(-1.2146981) q[2];
sx q[2];
rz(3.0991128) q[2];
rz(-2.5111607) q[3];
sx q[3];
rz(-2.5094331) q[3];
sx q[3];
rz(-2.650034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7522488) q[0];
sx q[0];
rz(-1.9135973) q[0];
sx q[0];
rz(0.4831627) q[0];
rz(2.0893611) q[1];
sx q[1];
rz(-1.1455043) q[1];
sx q[1];
rz(2.5767456) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42445499) q[0];
sx q[0];
rz(-1.0466482) q[0];
sx q[0];
rz(-2.9852887) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0658423) q[2];
sx q[2];
rz(-1.9905914) q[2];
sx q[2];
rz(2.9043353) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9996623) q[1];
sx q[1];
rz(-1.6301486) q[1];
sx q[1];
rz(1.1430986) q[1];
rz(-pi) q[2];
rz(2.5552093) q[3];
sx q[3];
rz(-0.27554232) q[3];
sx q[3];
rz(1.9540209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5715282) q[2];
sx q[2];
rz(-2.0740985) q[2];
sx q[2];
rz(1.8035536) q[2];
rz(1.3048874) q[3];
sx q[3];
rz(-2.0740502) q[3];
sx q[3];
rz(2.4664972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9682482) q[0];
sx q[0];
rz(-1.4302379) q[0];
sx q[0];
rz(0.54779732) q[0];
rz(-0.785218) q[1];
sx q[1];
rz(-1.806587) q[1];
sx q[1];
rz(-0.26842591) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9831532) q[0];
sx q[0];
rz(-1.325325) q[0];
sx q[0];
rz(-1.4863187) q[0];
x q[1];
rz(0.7873017) q[2];
sx q[2];
rz(-2.1848218) q[2];
sx q[2];
rz(1.8067443) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7859898) q[1];
sx q[1];
rz(-2.5463748) q[1];
sx q[1];
rz(1.3952414) q[1];
rz(-pi) q[2];
rz(2.1106342) q[3];
sx q[3];
rz(-1.3098048) q[3];
sx q[3];
rz(-0.48418448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.61775529) q[2];
sx q[2];
rz(-0.80344168) q[2];
sx q[2];
rz(-0.8141554) q[2];
rz(0.37627775) q[3];
sx q[3];
rz(-1.1637996) q[3];
sx q[3];
rz(0.072908727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5557264) q[0];
sx q[0];
rz(-1.4050452) q[0];
sx q[0];
rz(-1.0193753) q[0];
rz(0.85340071) q[1];
sx q[1];
rz(-1.1420206) q[1];
sx q[1];
rz(2.6928435) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9578581) q[0];
sx q[0];
rz(-1.0976037) q[0];
sx q[0];
rz(-2.8790022) q[0];
rz(2.8837187) q[2];
sx q[2];
rz(-1.2784625) q[2];
sx q[2];
rz(2.984798) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6905578) q[1];
sx q[1];
rz(-1.0447377) q[1];
sx q[1];
rz(1.2760389) q[1];
rz(-0.87047808) q[3];
sx q[3];
rz(-2.212489) q[3];
sx q[3];
rz(-0.93922797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2723508) q[2];
sx q[2];
rz(-1.3808455) q[2];
sx q[2];
rz(2.8273919) q[2];
rz(0.82434404) q[3];
sx q[3];
rz(-2.6894675) q[3];
sx q[3];
rz(0.79469386) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0712414) q[0];
sx q[0];
rz(-3.0817139) q[0];
sx q[0];
rz(-1.8810133) q[0];
rz(-2.4977327) q[1];
sx q[1];
rz(-1.9088129) q[1];
sx q[1];
rz(3.1226645) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8185794) q[0];
sx q[0];
rz(-1.5886663) q[0];
sx q[0];
rz(2.8779526) q[0];
rz(2.5817714) q[2];
sx q[2];
rz(-2.9580742) q[2];
sx q[2];
rz(-0.57932094) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.85266528) q[1];
sx q[1];
rz(-1.7551433) q[1];
sx q[1];
rz(1.5194555) q[1];
rz(-1.3350248) q[3];
sx q[3];
rz(-1.2887495) q[3];
sx q[3];
rz(2.3896133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5921322) q[2];
sx q[2];
rz(-0.38828725) q[2];
sx q[2];
rz(-0.67031676) q[2];
rz(2.629225) q[3];
sx q[3];
rz(-1.3918326) q[3];
sx q[3];
rz(1.4204773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(0.55384127) q[0];
sx q[0];
rz(-1.8974263) q[0];
sx q[0];
rz(0.49945369) q[0];
rz(1.5746501) q[1];
sx q[1];
rz(-2.8630239) q[1];
sx q[1];
rz(-2.0589028) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2912746) q[0];
sx q[0];
rz(-0.58915888) q[0];
sx q[0];
rz(-0.088081443) q[0];
rz(-pi) q[1];
rz(-2.322299) q[2];
sx q[2];
rz(-1.5466006) q[2];
sx q[2];
rz(2.1404612) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.06304) q[1];
sx q[1];
rz(-1.5895956) q[1];
sx q[1];
rz(0.54363721) q[1];
x q[2];
rz(1.2438251) q[3];
sx q[3];
rz(-1.9332814) q[3];
sx q[3];
rz(2.8231951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3988951) q[2];
sx q[2];
rz(-1.9766786) q[2];
sx q[2];
rz(2.6297074) q[2];
rz(0.39294696) q[3];
sx q[3];
rz(-1.4018551) q[3];
sx q[3];
rz(-2.0846562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.186541) q[0];
sx q[0];
rz(-1.3165836) q[0];
sx q[0];
rz(-2.5008428) q[0];
rz(2.4004249) q[1];
sx q[1];
rz(-2.3186431) q[1];
sx q[1];
rz(2.9021312) q[1];
rz(-2.1869833) q[2];
sx q[2];
rz(-1.2381427) q[2];
sx q[2];
rz(2.8340813) q[2];
rz(-2.3951204) q[3];
sx q[3];
rz(-1.4440047) q[3];
sx q[3];
rz(0.11228893) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
