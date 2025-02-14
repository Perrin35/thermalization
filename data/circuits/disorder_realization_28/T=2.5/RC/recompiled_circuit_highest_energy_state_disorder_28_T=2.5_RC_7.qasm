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
rz(-2.3856491) q[0];
sx q[0];
rz(-0.33060253) q[0];
sx q[0];
rz(0.27275738) q[0];
rz(2.7859712) q[1];
sx q[1];
rz(-2.1115117) q[1];
sx q[1];
rz(-1.570809) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8518675) q[0];
sx q[0];
rz(-0.12061943) q[0];
sx q[0];
rz(0.9319181) q[0];
rz(-pi) q[1];
rz(-2.5433538) q[2];
sx q[2];
rz(-0.5454059) q[2];
sx q[2];
rz(2.0871833) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.570191) q[1];
sx q[1];
rz(-0.68059151) q[1];
sx q[1];
rz(0.067318214) q[1];
rz(-1.4727888) q[3];
sx q[3];
rz(-1.2036754) q[3];
sx q[3];
rz(-0.29686061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1640799) q[2];
sx q[2];
rz(-1.7089607) q[2];
sx q[2];
rz(-2.254503) q[2];
rz(1.2648434) q[3];
sx q[3];
rz(-1.0586459) q[3];
sx q[3];
rz(-3.1177055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.7510659) q[0];
sx q[0];
rz(-2.2845415) q[0];
sx q[0];
rz(-2.8476025) q[0];
rz(-1.7901621) q[1];
sx q[1];
rz(-1.0963115) q[1];
sx q[1];
rz(1.1658123) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4664) q[0];
sx q[0];
rz(-1.8972016) q[0];
sx q[0];
rz(-2.9221888) q[0];
rz(-1.2109257) q[2];
sx q[2];
rz(-3.0635298) q[2];
sx q[2];
rz(0.087457267) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9820751) q[1];
sx q[1];
rz(-1.0613835) q[1];
sx q[1];
rz(1.0647091) q[1];
rz(-0.43710093) q[3];
sx q[3];
rz(-0.61449285) q[3];
sx q[3];
rz(-1.3853427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4531648) q[2];
sx q[2];
rz(-2.5025949) q[2];
sx q[2];
rz(-1.9387908) q[2];
rz(-1.5419675) q[3];
sx q[3];
rz(-1.8034233) q[3];
sx q[3];
rz(-0.28551027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0315392) q[0];
sx q[0];
rz(-3.0634614) q[0];
sx q[0];
rz(-1.9670638) q[0];
rz(0.9749167) q[1];
sx q[1];
rz(-2.2623623) q[1];
sx q[1];
rz(2.6851795) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83041265) q[0];
sx q[0];
rz(-1.6025271) q[0];
sx q[0];
rz(-1.8422442) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1295616) q[2];
sx q[2];
rz(-0.44545275) q[2];
sx q[2];
rz(0.73778668) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.76589903) q[1];
sx q[1];
rz(-1.0119034) q[1];
sx q[1];
rz(0.8984579) q[1];
rz(-pi) q[2];
rz(2.1022845) q[3];
sx q[3];
rz(-0.78954299) q[3];
sx q[3];
rz(1.2293775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5259033) q[2];
sx q[2];
rz(-2.758785) q[2];
sx q[2];
rz(2.3846386) q[2];
rz(3.1225539) q[3];
sx q[3];
rz(-1.8502539) q[3];
sx q[3];
rz(-0.835787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2794936) q[0];
sx q[0];
rz(-2.4308496) q[0];
sx q[0];
rz(-0.94974649) q[0];
rz(-0.21385916) q[1];
sx q[1];
rz(-2.2016134) q[1];
sx q[1];
rz(0.59923879) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.43219) q[0];
sx q[0];
rz(-2.4718542) q[0];
sx q[0];
rz(1.7959005) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0867351) q[2];
sx q[2];
rz(-1.4241059) q[2];
sx q[2];
rz(-0.88146081) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9436151) q[1];
sx q[1];
rz(-0.57398263) q[1];
sx q[1];
rz(1.9616433) q[1];
rz(1.7091126) q[3];
sx q[3];
rz(-0.99405046) q[3];
sx q[3];
rz(0.53046295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5102101) q[2];
sx q[2];
rz(-1.5082794) q[2];
sx q[2];
rz(-2.317826) q[2];
rz(-1.0907762) q[3];
sx q[3];
rz(-2.3407276) q[3];
sx q[3];
rz(2.4763079) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7432231) q[0];
sx q[0];
rz(-2.7581765) q[0];
sx q[0];
rz(-0.98633352) q[0];
rz(-0.24364722) q[1];
sx q[1];
rz(-1.8758834) q[1];
sx q[1];
rz(2.9820014) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0054609912) q[0];
sx q[0];
rz(-1.8842995) q[0];
sx q[0];
rz(0.82833293) q[0];
x q[1];
rz(-0.56964376) q[2];
sx q[2];
rz(-2.2692371) q[2];
sx q[2];
rz(-1.6676355) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1732542) q[1];
sx q[1];
rz(-2.6880777) q[1];
sx q[1];
rz(-1.7733002) q[1];
rz(0.29496622) q[3];
sx q[3];
rz(-2.2497503) q[3];
sx q[3];
rz(-2.7243638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7470982) q[2];
sx q[2];
rz(-1.8793841) q[2];
sx q[2];
rz(-2.0060284) q[2];
rz(0.46843946) q[3];
sx q[3];
rz(-1.2218385) q[3];
sx q[3];
rz(0.92456094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0202476) q[0];
sx q[0];
rz(-2.0609042) q[0];
sx q[0];
rz(0.24170804) q[0];
rz(-1.3555591) q[1];
sx q[1];
rz(-2.2691998) q[1];
sx q[1];
rz(0.88368574) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51118219) q[0];
sx q[0];
rz(-0.97994655) q[0];
sx q[0];
rz(0.972959) q[0];
rz(-1.3140129) q[2];
sx q[2];
rz(-1.2681539) q[2];
sx q[2];
rz(0.0079449991) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.984772) q[1];
sx q[1];
rz(-1.7659016) q[1];
sx q[1];
rz(-2.0181993) q[1];
rz(-pi) q[2];
rz(-3.1413636) q[3];
sx q[3];
rz(-1.9817686) q[3];
sx q[3];
rz(0.64490333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9200865) q[2];
sx q[2];
rz(-0.93044996) q[2];
sx q[2];
rz(2.7152756) q[2];
rz(0.91655556) q[3];
sx q[3];
rz(-1.5896268) q[3];
sx q[3];
rz(-1.1030654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78974462) q[0];
sx q[0];
rz(-1.177657) q[0];
sx q[0];
rz(2.015693) q[0];
rz(-0.060983505) q[1];
sx q[1];
rz(-0.95044249) q[1];
sx q[1];
rz(1.0152063) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7495817) q[0];
sx q[0];
rz(-0.77031743) q[0];
sx q[0];
rz(2.5213741) q[0];
x q[1];
rz(-0.073748579) q[2];
sx q[2];
rz(-1.2332067) q[2];
sx q[2];
rz(0.13843564) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0132844) q[1];
sx q[1];
rz(-0.40402141) q[1];
sx q[1];
rz(-1.4133417) q[1];
x q[2];
rz(3.0281248) q[3];
sx q[3];
rz(-1.5312636) q[3];
sx q[3];
rz(3.0276445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.31275493) q[2];
sx q[2];
rz(-1.2042896) q[2];
sx q[2];
rz(0.22605669) q[2];
rz(-0.95212805) q[3];
sx q[3];
rz(-0.13855562) q[3];
sx q[3];
rz(-0.8086732) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011618135) q[0];
sx q[0];
rz(-1.3148774) q[0];
sx q[0];
rz(-2.8011978) q[0];
rz(-3.0420692) q[1];
sx q[1];
rz(-2.0390022) q[1];
sx q[1];
rz(1.5460825) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81170245) q[0];
sx q[0];
rz(-1.2445589) q[0];
sx q[0];
rz(0.6321815) q[0];
rz(1.298102) q[2];
sx q[2];
rz(-1.6249949) q[2];
sx q[2];
rz(1.3860089) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.73879646) q[1];
sx q[1];
rz(-1.6127043) q[1];
sx q[1];
rz(-0.099385029) q[1];
x q[2];
rz(-2.9603954) q[3];
sx q[3];
rz(-2.6799723) q[3];
sx q[3];
rz(0.68440729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.16326363) q[2];
sx q[2];
rz(-1.1922057) q[2];
sx q[2];
rz(1.0015798) q[2];
rz(-1.7892276) q[3];
sx q[3];
rz(-2.6632301) q[3];
sx q[3];
rz(1.1866239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1773862) q[0];
sx q[0];
rz(-2.0851676) q[0];
sx q[0];
rz(0.00038432234) q[0];
rz(0.43801019) q[1];
sx q[1];
rz(-1.6338467) q[1];
sx q[1];
rz(-0.45723525) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9033636) q[0];
sx q[0];
rz(-1.6203124) q[0];
sx q[0];
rz(-1.7286506) q[0];
x q[1];
rz(-2.5633144) q[2];
sx q[2];
rz(-2.1077029) q[2];
sx q[2];
rz(1.121305) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.286736) q[1];
sx q[1];
rz(-1.348362) q[1];
sx q[1];
rz(-0.35448928) q[1];
x q[2];
rz(-0.5064241) q[3];
sx q[3];
rz(-1.862414) q[3];
sx q[3];
rz(-2.9675561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.75927258) q[2];
sx q[2];
rz(-0.98733941) q[2];
sx q[2];
rz(-1.5116723) q[2];
rz(3.0509389) q[3];
sx q[3];
rz(-2.1001308) q[3];
sx q[3];
rz(-0.67971984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35336211) q[0];
sx q[0];
rz(-1.7857977) q[0];
sx q[0];
rz(-0.28025383) q[0];
rz(1.7578112) q[1];
sx q[1];
rz(-2.6774112) q[1];
sx q[1];
rz(-1.9162477) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2458151) q[0];
sx q[0];
rz(-1.7208463) q[0];
sx q[0];
rz(1.4729985) q[0];
rz(-pi) q[1];
rz(2.2544075) q[2];
sx q[2];
rz(-1.9211743) q[2];
sx q[2];
rz(-0.35752192) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.30825057) q[1];
sx q[1];
rz(-2.3099515) q[1];
sx q[1];
rz(1.2396481) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8605385) q[3];
sx q[3];
rz(-1.8612766) q[3];
sx q[3];
rz(-1.2742467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9489991) q[2];
sx q[2];
rz(-1.6152363) q[2];
sx q[2];
rz(-2.4601649) q[2];
rz(-1.1000018) q[3];
sx q[3];
rz(-2.2097578) q[3];
sx q[3];
rz(-1.8345691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6901907) q[0];
sx q[0];
rz(-0.9342397) q[0];
sx q[0];
rz(-1.1697212) q[0];
rz(-2.7611217) q[1];
sx q[1];
rz(-2.4741551) q[1];
sx q[1];
rz(-2.9173775) q[1];
rz(-1.8809594) q[2];
sx q[2];
rz(-0.72523965) q[2];
sx q[2];
rz(0.88015866) q[2];
rz(-2.1721645) q[3];
sx q[3];
rz(-2.457544) q[3];
sx q[3];
rz(-2.8586879) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
