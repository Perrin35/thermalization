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
rz(2.9444115) q[0];
sx q[0];
rz(-1.0721167) q[0];
sx q[0];
rz(-1.5972692) q[0];
rz(0.3577258) q[1];
sx q[1];
rz(-1.3074713) q[1];
sx q[1];
rz(0.37582418) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2471814) q[0];
sx q[0];
rz(-0.80601776) q[0];
sx q[0];
rz(-1.4708118) q[0];
rz(1.8456468) q[2];
sx q[2];
rz(-1.0604727) q[2];
sx q[2];
rz(-2.9269671) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.5511939) q[1];
sx q[1];
rz(-2.2839468) q[1];
sx q[1];
rz(0.80618315) q[1];
rz(1.1268257) q[3];
sx q[3];
rz(-0.68279632) q[3];
sx q[3];
rz(-1.7449774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4481005) q[2];
sx q[2];
rz(-0.38303146) q[2];
sx q[2];
rz(-2.7310272) q[2];
rz(0.52452123) q[3];
sx q[3];
rz(-0.21206847) q[3];
sx q[3];
rz(-1.9922403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18325026) q[0];
sx q[0];
rz(-0.14160937) q[0];
sx q[0];
rz(-0.70567065) q[0];
rz(-1.2880098) q[1];
sx q[1];
rz(-2.5649004) q[1];
sx q[1];
rz(-1.1340595) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9767446) q[0];
sx q[0];
rz(-2.3253092) q[0];
sx q[0];
rz(-1.0036941) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7984624) q[2];
sx q[2];
rz(-1.742082) q[2];
sx q[2];
rz(-0.73395455) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6780939) q[1];
sx q[1];
rz(-1.948852) q[1];
sx q[1];
rz(2.3854756) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5763578) q[3];
sx q[3];
rz(-1.4197403) q[3];
sx q[3];
rz(-0.24542576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3386393) q[2];
sx q[2];
rz(-2.6889375) q[2];
sx q[2];
rz(0.53042859) q[2];
rz(0.28702921) q[3];
sx q[3];
rz(-0.48873264) q[3];
sx q[3];
rz(0.73634017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1257783) q[0];
sx q[0];
rz(-2.0398572) q[0];
sx q[0];
rz(-1.0544448) q[0];
rz(-1.6619445) q[1];
sx q[1];
rz(-0.92250657) q[1];
sx q[1];
rz(1.066801) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22196348) q[0];
sx q[0];
rz(-2.9502075) q[0];
sx q[0];
rz(-1.3758977) q[0];
x q[1];
rz(-2.4652723) q[2];
sx q[2];
rz(-1.6427077) q[2];
sx q[2];
rz(-1.423171) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.57745614) q[1];
sx q[1];
rz(-1.3423663) q[1];
sx q[1];
rz(-1.0962568) q[1];
rz(-pi) q[2];
rz(2.2852632) q[3];
sx q[3];
rz(-1.3420055) q[3];
sx q[3];
rz(-1.8280054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.58831424) q[2];
sx q[2];
rz(-1.8411627) q[2];
sx q[2];
rz(1.1021357) q[2];
rz(-0.45840248) q[3];
sx q[3];
rz(-1.6130684) q[3];
sx q[3];
rz(2.5100191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.691953) q[0];
sx q[0];
rz(-0.75943685) q[0];
sx q[0];
rz(2.6547883) q[0];
rz(1.5961157) q[1];
sx q[1];
rz(-0.38914746) q[1];
sx q[1];
rz(2.5130689) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8299943) q[0];
sx q[0];
rz(-1.4382558) q[0];
sx q[0];
rz(1.1740985) q[0];
x q[1];
rz(-0.14735584) q[2];
sx q[2];
rz(-1.0026081) q[2];
sx q[2];
rz(-2.1939366) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6316136) q[1];
sx q[1];
rz(-1.1812796) q[1];
sx q[1];
rz(-2.9385467) q[1];
x q[2];
rz(-1.6877141) q[3];
sx q[3];
rz(-0.55171574) q[3];
sx q[3];
rz(-1.6601613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.57807606) q[2];
sx q[2];
rz(-3.0264644) q[2];
sx q[2];
rz(0.76187491) q[2];
rz(0.21841194) q[3];
sx q[3];
rz(-0.29774791) q[3];
sx q[3];
rz(2.0930321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31576425) q[0];
sx q[0];
rz(-2.789412) q[0];
sx q[0];
rz(2.5370362) q[0];
rz(-1.293659) q[1];
sx q[1];
rz(-0.81984055) q[1];
sx q[1];
rz(-2.7972319) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47681706) q[0];
sx q[0];
rz(-1.6533543) q[0];
sx q[0];
rz(0.020391012) q[0];
rz(1.2372962) q[2];
sx q[2];
rz(-0.97003675) q[2];
sx q[2];
rz(-1.5983901) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1188914) q[1];
sx q[1];
rz(-1.1384607) q[1];
sx q[1];
rz(-2.6328264) q[1];
rz(-pi) q[2];
rz(0.46724288) q[3];
sx q[3];
rz(-2.1911216) q[3];
sx q[3];
rz(1.0091227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0209501) q[2];
sx q[2];
rz(-0.63434333) q[2];
sx q[2];
rz(-2.9316588) q[2];
rz(1.1076934) q[3];
sx q[3];
rz(-1.6898797) q[3];
sx q[3];
rz(-0.2618739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8629859) q[0];
sx q[0];
rz(-0.59026533) q[0];
sx q[0];
rz(3.0721831) q[0];
rz(-0.033992652) q[1];
sx q[1];
rz(-2.4885204) q[1];
sx q[1];
rz(0.45190826) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6959608) q[0];
sx q[0];
rz(-1.9602141) q[0];
sx q[0];
rz(-1.7444064) q[0];
rz(0.71058013) q[2];
sx q[2];
rz(-2.3167774) q[2];
sx q[2];
rz(-1.9503066) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.76594162) q[1];
sx q[1];
rz(-0.022863511) q[1];
sx q[1];
rz(1.8909062) q[1];
rz(-2.8699401) q[3];
sx q[3];
rz(-1.8166887) q[3];
sx q[3];
rz(0.065448381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7273093) q[2];
sx q[2];
rz(-2.3392129) q[2];
sx q[2];
rz(-1.9963973) q[2];
rz(0.9582054) q[3];
sx q[3];
rz(-1.9327791) q[3];
sx q[3];
rz(-2.8632979) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77022839) q[0];
sx q[0];
rz(-0.48786491) q[0];
sx q[0];
rz(2.2801953) q[0];
rz(-0.96249145) q[1];
sx q[1];
rz(-0.95314127) q[1];
sx q[1];
rz(2.9846587) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3478137) q[0];
sx q[0];
rz(-2.2484178) q[0];
sx q[0];
rz(0.47394808) q[0];
rz(-pi) q[1];
rz(-0.10617039) q[2];
sx q[2];
rz(-1.9261126) q[2];
sx q[2];
rz(0.67007845) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7879856) q[1];
sx q[1];
rz(-0.79670409) q[1];
sx q[1];
rz(0.70938797) q[1];
rz(-pi) q[2];
rz(1.3577363) q[3];
sx q[3];
rz(-1.5231265) q[3];
sx q[3];
rz(-0.21040645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.12338766) q[2];
sx q[2];
rz(-2.3473479) q[2];
sx q[2];
rz(2.0920853) q[2];
rz(-2.6617995) q[3];
sx q[3];
rz(-2.9496461) q[3];
sx q[3];
rz(-0.16201365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70742339) q[0];
sx q[0];
rz(-0.24288414) q[0];
sx q[0];
rz(2.6450787) q[0];
rz(1.4855344) q[1];
sx q[1];
rz(-0.86692202) q[1];
sx q[1];
rz(3.1190125) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17538769) q[0];
sx q[0];
rz(-1.3577355) q[0];
sx q[0];
rz(-2.9609462) q[0];
rz(-pi) q[1];
rz(1.1316001) q[2];
sx q[2];
rz(-0.15838693) q[2];
sx q[2];
rz(1.8902376) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.034790967) q[1];
sx q[1];
rz(-0.80134922) q[1];
sx q[1];
rz(-1.31557) q[1];
rz(1.0362421) q[3];
sx q[3];
rz(-1.4009985) q[3];
sx q[3];
rz(-1.9056729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.61985832) q[2];
sx q[2];
rz(-1.1405742) q[2];
sx q[2];
rz(0.015259585) q[2];
rz(-0.36733019) q[3];
sx q[3];
rz(-0.44706774) q[3];
sx q[3];
rz(-2.7801311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.5527282) q[0];
sx q[0];
rz(-0.22607729) q[0];
sx q[0];
rz(0.077762522) q[0];
rz(3.0231158) q[1];
sx q[1];
rz(-2.8006554) q[1];
sx q[1];
rz(-0.65611398) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2099119) q[0];
sx q[0];
rz(-2.8114578) q[0];
sx q[0];
rz(-1.1709377) q[0];
rz(-0.17123789) q[2];
sx q[2];
rz(-1.3067553) q[2];
sx q[2];
rz(-1.1080826) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.205913) q[1];
sx q[1];
rz(-0.64710837) q[1];
sx q[1];
rz(0.86828219) q[1];
rz(-pi) q[2];
rz(-1.0540851) q[3];
sx q[3];
rz(-1.0789405) q[3];
sx q[3];
rz(-2.3824177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6180827) q[2];
sx q[2];
rz(-1.9205576) q[2];
sx q[2];
rz(2.6005884) q[2];
rz(-2.6909761) q[3];
sx q[3];
rz(-1.6559947) q[3];
sx q[3];
rz(0.97644067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70380521) q[0];
sx q[0];
rz(-1.7923651) q[0];
sx q[0];
rz(-0.10677862) q[0];
rz(1.59185) q[1];
sx q[1];
rz(-0.50250643) q[1];
sx q[1];
rz(1.3776616) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6984682) q[0];
sx q[0];
rz(-1.6148991) q[0];
sx q[0];
rz(-1.1209784) q[0];
rz(0.3129199) q[2];
sx q[2];
rz(-1.358479) q[2];
sx q[2];
rz(0.91845271) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.622887) q[1];
sx q[1];
rz(-2.7939502) q[1];
sx q[1];
rz(-0.80153377) q[1];
x q[2];
rz(-0.53394188) q[3];
sx q[3];
rz(-2.0755872) q[3];
sx q[3];
rz(1.8478888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0489073) q[2];
sx q[2];
rz(-1.4212298) q[2];
sx q[2];
rz(0.29472026) q[2];
rz(-2.5367391) q[3];
sx q[3];
rz(-1.2302783) q[3];
sx q[3];
rz(-3.0321583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7034364) q[0];
sx q[0];
rz(-1.823517) q[0];
sx q[0];
rz(2.1448268) q[0];
rz(-3.1238212) q[1];
sx q[1];
rz(-1.8780864) q[1];
sx q[1];
rz(-1.3395739) q[1];
rz(2.3370635) q[2];
sx q[2];
rz(-1.1297516) q[2];
sx q[2];
rz(-0.17937112) q[2];
rz(-2.3066219) q[3];
sx q[3];
rz(-2.2899262) q[3];
sx q[3];
rz(-0.038577608) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
