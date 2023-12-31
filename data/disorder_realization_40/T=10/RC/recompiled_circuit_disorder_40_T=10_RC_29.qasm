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
rz(-1.4594266) q[1];
sx q[1];
rz(-1.6571801) q[1];
sx q[1];
rz(0.15375528) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0271887) q[0];
sx q[0];
rz(-1.1429042) q[0];
sx q[0];
rz(1.2501636) q[0];
rz(-pi) q[1];
rz(-1.8122187) q[2];
sx q[2];
rz(-1.6454988) q[2];
sx q[2];
rz(2.9818997) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1195214) q[1];
sx q[1];
rz(-1.6011366) q[1];
sx q[1];
rz(-1.2519217) q[1];
rz(2.0102242) q[3];
sx q[3];
rz(-1.8779552) q[3];
sx q[3];
rz(0.22104056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.443632) q[2];
sx q[2];
rz(-1.4322832) q[2];
sx q[2];
rz(1.704818) q[2];
rz(0.73389655) q[3];
sx q[3];
rz(-1.5489483) q[3];
sx q[3];
rz(-2.6255887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2519418) q[0];
sx q[0];
rz(-1.2263068) q[0];
sx q[0];
rz(-0.92457986) q[0];
rz(0.997116) q[1];
sx q[1];
rz(-2.6328502) q[1];
sx q[1];
rz(1.8181713) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88432089) q[0];
sx q[0];
rz(-1.1182251) q[0];
sx q[0];
rz(1.0898468) q[0];
rz(-0.3496062) q[2];
sx q[2];
rz(-1.5290302) q[2];
sx q[2];
rz(-1.8984399) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.70222774) q[1];
sx q[1];
rz(-0.80184466) q[1];
sx q[1];
rz(2.761809) q[1];
x q[2];
rz(-0.38815659) q[3];
sx q[3];
rz(-1.0309439) q[3];
sx q[3];
rz(-0.29263465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.20415846) q[2];
sx q[2];
rz(-1.5896475) q[2];
sx q[2];
rz(2.3834174) q[2];
rz(-2.5126863) q[3];
sx q[3];
rz(-2.7401676) q[3];
sx q[3];
rz(-1.988407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73308289) q[0];
sx q[0];
rz(-1.2671616) q[0];
sx q[0];
rz(2.4531903) q[0];
rz(3.0738661) q[1];
sx q[1];
rz(-1.3893145) q[1];
sx q[1];
rz(0.53007954) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.006223) q[0];
sx q[0];
rz(-2.831499) q[0];
sx q[0];
rz(-1.2139411) q[0];
rz(-pi) q[1];
rz(-2.0968139) q[2];
sx q[2];
rz(-1.194343) q[2];
sx q[2];
rz(-1.2029755) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9024132) q[1];
sx q[1];
rz(-1.2352408) q[1];
sx q[1];
rz(-2.2526342) q[1];
rz(-pi) q[2];
rz(-1.1881371) q[3];
sx q[3];
rz(-0.2573765) q[3];
sx q[3];
rz(2.5352258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3337341) q[2];
sx q[2];
rz(-0.01161751) q[2];
sx q[2];
rz(-0.90144908) q[2];
rz(0.83550134) q[3];
sx q[3];
rz(-1.52799) q[3];
sx q[3];
rz(-1.3114595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19105844) q[0];
sx q[0];
rz(-1.598851) q[0];
sx q[0];
rz(2.3572671) q[0];
rz(-0.061231881) q[1];
sx q[1];
rz(-0.71413723) q[1];
sx q[1];
rz(3.004946) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.559583) q[0];
sx q[0];
rz(-2.0499381) q[0];
sx q[0];
rz(-2.9942306) q[0];
rz(-pi) q[1];
rz(-2.5606974) q[2];
sx q[2];
rz(-2.6022606) q[2];
sx q[2];
rz(0.5667333) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5977051) q[1];
sx q[1];
rz(-1.502431) q[1];
sx q[1];
rz(-3.1176561) q[1];
rz(-pi) q[2];
rz(3.0487719) q[3];
sx q[3];
rz(-2.1448359) q[3];
sx q[3];
rz(1.012158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1293929) q[2];
sx q[2];
rz(-0.94038525) q[2];
sx q[2];
rz(-2.5811035) q[2];
rz(-3.1292606) q[3];
sx q[3];
rz(-2.2380232) q[3];
sx q[3];
rz(1.0906609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.7085768) q[0];
sx q[0];
rz(-2.5362159) q[0];
sx q[0];
rz(0.82114712) q[0];
rz(-0.87617809) q[1];
sx q[1];
rz(-0.89996243) q[1];
sx q[1];
rz(-1.7339773) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9579904) q[0];
sx q[0];
rz(-1.9200268) q[0];
sx q[0];
rz(2.7244085) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1528835) q[2];
sx q[2];
rz(-1.313813) q[2];
sx q[2];
rz(2.7866521) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.69406063) q[1];
sx q[1];
rz(-2.5563055) q[1];
sx q[1];
rz(1.3259757) q[1];
rz(-pi) q[2];
rz(-1.9305265) q[3];
sx q[3];
rz(-1.3621646) q[3];
sx q[3];
rz(-1.2071351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.451482) q[2];
sx q[2];
rz(-1.9268945) q[2];
sx q[2];
rz(3.0991128) q[2];
rz(-2.5111607) q[3];
sx q[3];
rz(-0.63215956) q[3];
sx q[3];
rz(2.650034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7522488) q[0];
sx q[0];
rz(-1.9135973) q[0];
sx q[0];
rz(-2.65843) q[0];
rz(-1.0522316) q[1];
sx q[1];
rz(-1.1455043) q[1];
sx q[1];
rz(-0.56484708) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4120867) q[0];
sx q[0];
rz(-2.5967263) q[0];
sx q[0];
rz(-1.3077523) q[0];
rz(-pi) q[1];
rz(1.9916612) q[2];
sx q[2];
rz(-1.5016342) q[2];
sx q[2];
rz(1.3644621) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40183345) q[1];
sx q[1];
rz(-1.1439011) q[1];
sx q[1];
rz(0.065211936) q[1];
x q[2];
rz(1.7259898) q[3];
sx q[3];
rz(-1.7994013) q[3];
sx q[3];
rz(0.58333635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5715282) q[2];
sx q[2];
rz(-1.0674942) q[2];
sx q[2];
rz(1.338039) q[2];
rz(-1.3048874) q[3];
sx q[3];
rz(-2.0740502) q[3];
sx q[3];
rz(-2.4664972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17334443) q[0];
sx q[0];
rz(-1.7113547) q[0];
sx q[0];
rz(2.5937953) q[0];
rz(-2.3563747) q[1];
sx q[1];
rz(-1.3350057) q[1];
sx q[1];
rz(2.8731667) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4329322) q[0];
sx q[0];
rz(-1.6527358) q[0];
sx q[0];
rz(2.8952778) q[0];
rz(-2.3586876) q[2];
sx q[2];
rz(-0.9559388) q[2];
sx q[2];
rz(0.2851141) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5749579) q[1];
sx q[1];
rz(-2.1556427) q[1];
sx q[1];
rz(0.11771867) q[1];
rz(-pi) q[2];
rz(-1.0309585) q[3];
sx q[3];
rz(-1.3098048) q[3];
sx q[3];
rz(2.6574082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.61775529) q[2];
sx q[2];
rz(-0.80344168) q[2];
sx q[2];
rz(0.8141554) q[2];
rz(-2.7653149) q[3];
sx q[3];
rz(-1.9777931) q[3];
sx q[3];
rz(-0.072908727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58586621) q[0];
sx q[0];
rz(-1.4050452) q[0];
sx q[0];
rz(1.0193753) q[0];
rz(0.85340071) q[1];
sx q[1];
rz(-1.9995721) q[1];
sx q[1];
rz(-2.6928435) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4067626) q[0];
sx q[0];
rz(-1.3376298) q[0];
sx q[0];
rz(-2.0582817) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8837187) q[2];
sx q[2];
rz(-1.2784625) q[2];
sx q[2];
rz(-0.1567947) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6905578) q[1];
sx q[1];
rz(-1.0447377) q[1];
sx q[1];
rz(1.2760389) q[1];
rz(-pi) q[2];
rz(2.2711146) q[3];
sx q[3];
rz(-2.212489) q[3];
sx q[3];
rz(-0.93922797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2723508) q[2];
sx q[2];
rz(-1.3808455) q[2];
sx q[2];
rz(-0.31420079) q[2];
rz(0.82434404) q[3];
sx q[3];
rz(-2.6894675) q[3];
sx q[3];
rz(-2.3468988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0712414) q[0];
sx q[0];
rz(-3.0817139) q[0];
sx q[0];
rz(1.8810133) q[0];
rz(0.64385995) q[1];
sx q[1];
rz(-1.9088129) q[1];
sx q[1];
rz(3.1226645) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8185794) q[0];
sx q[0];
rz(-1.5886663) q[0];
sx q[0];
rz(0.26364003) q[0];
rz(1.47255) q[2];
sx q[2];
rz(-1.7260523) q[2];
sx q[2];
rz(1.1467903) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70871204) q[1];
sx q[1];
rz(-1.5203262) q[1];
sx q[1];
rz(0.18458472) q[1];
x q[2];
rz(0.67846672) q[3];
sx q[3];
rz(-2.7760091) q[3];
sx q[3];
rz(-1.4640704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.54946047) q[2];
sx q[2];
rz(-0.38828725) q[2];
sx q[2];
rz(2.4712759) q[2];
rz(2.629225) q[3];
sx q[3];
rz(-1.3918326) q[3];
sx q[3];
rz(-1.7211154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55384127) q[0];
sx q[0];
rz(-1.2441664) q[0];
sx q[0];
rz(0.49945369) q[0];
rz(-1.5669426) q[1];
sx q[1];
rz(-0.27856871) q[1];
sx q[1];
rz(2.0589028) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1854541) q[0];
sx q[0];
rz(-2.1573665) q[0];
sx q[0];
rz(-1.5120718) q[0];
rz(-3.1084836) q[2];
sx q[2];
rz(-2.3220255) q[2];
sx q[2];
rz(-0.59226743) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4808828) q[1];
sx q[1];
rz(-2.1143267) q[1];
sx q[1];
rz(-1.5927614) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8977676) q[3];
sx q[3];
rz(-1.9332814) q[3];
sx q[3];
rz(2.8231951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7426976) q[2];
sx q[2];
rz(-1.9766786) q[2];
sx q[2];
rz(0.51188525) q[2];
rz(0.39294696) q[3];
sx q[3];
rz(-1.4018551) q[3];
sx q[3];
rz(-2.0846562) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.186541) q[0];
sx q[0];
rz(-1.3165836) q[0];
sx q[0];
rz(-2.5008428) q[0];
rz(-0.74116771) q[1];
sx q[1];
rz(-2.3186431) q[1];
sx q[1];
rz(2.9021312) q[1];
rz(-1.0319866) q[2];
sx q[2];
rz(-2.451755) q[2];
sx q[2];
rz(-2.3103726) q[2];
rz(2.3951204) q[3];
sx q[3];
rz(-1.6975879) q[3];
sx q[3];
rz(-3.0293037) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
