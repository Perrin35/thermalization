OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.32970348) q[0];
sx q[0];
rz(-2.831037) q[0];
sx q[0];
rz(-1.1986873) q[0];
rz(-2.0074453) q[1];
sx q[1];
rz(2.3448047) q[1];
sx q[1];
rz(12.264836) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5422259) q[0];
sx q[0];
rz(-0.24655534) q[0];
sx q[0];
rz(-0.85041468) q[0];
x q[1];
rz(-1.0204591) q[2];
sx q[2];
rz(-2.7937903) q[2];
sx q[2];
rz(-0.90711275) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.63362279) q[1];
sx q[1];
rz(-2.3862641) q[1];
sx q[1];
rz(-1.669674) q[1];
rz(1.1687247) q[3];
sx q[3];
rz(-1.744606) q[3];
sx q[3];
rz(0.92648348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.6887168) q[2];
sx q[2];
rz(-1.0998925) q[2];
sx q[2];
rz(2.0067046) q[2];
rz(-3.0196043) q[3];
sx q[3];
rz(-2.205409) q[3];
sx q[3];
rz(0.056099135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4695456) q[0];
sx q[0];
rz(-2.8944954) q[0];
sx q[0];
rz(0.0023512996) q[0];
rz(-0.40070847) q[1];
sx q[1];
rz(-1.7727163) q[1];
sx q[1];
rz(-0.8561264) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24361783) q[0];
sx q[0];
rz(-1.8254733) q[0];
sx q[0];
rz(-0.93230482) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0547423) q[2];
sx q[2];
rz(-1.8212089) q[2];
sx q[2];
rz(2.5977787) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.90710956) q[1];
sx q[1];
rz(-2.6119348) q[1];
sx q[1];
rz(1.6698599) q[1];
rz(-0.54243954) q[3];
sx q[3];
rz(-2.503201) q[3];
sx q[3];
rz(0.29104376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1043642) q[2];
sx q[2];
rz(-1.9375485) q[2];
sx q[2];
rz(-0.47002235) q[2];
rz(0.59703279) q[3];
sx q[3];
rz(-3.0266914) q[3];
sx q[3];
rz(-1.443583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7715348) q[0];
sx q[0];
rz(-2.6486588) q[0];
sx q[0];
rz(-0.22802995) q[0];
rz(2.6920964) q[1];
sx q[1];
rz(-1.0003961) q[1];
sx q[1];
rz(2.1953348) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11433521) q[0];
sx q[0];
rz(-1.1722306) q[0];
sx q[0];
rz(-0.27862676) q[0];
x q[1];
rz(-2.3870638) q[2];
sx q[2];
rz(-2.3268394) q[2];
sx q[2];
rz(2.192564) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0277532) q[1];
sx q[1];
rz(-1.3885173) q[1];
sx q[1];
rz(-0.34232847) q[1];
rz(-pi) q[2];
rz(-2.0722996) q[3];
sx q[3];
rz(-2.8208102) q[3];
sx q[3];
rz(2.4298885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9049282) q[2];
sx q[2];
rz(-1.4651352) q[2];
sx q[2];
rz(1.6620592) q[2];
rz(-2.9605401) q[3];
sx q[3];
rz(-0.79020399) q[3];
sx q[3];
rz(-0.886206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7774696) q[0];
sx q[0];
rz(-1.5622219) q[0];
sx q[0];
rz(0.89853483) q[0];
rz(-2.6354375) q[1];
sx q[1];
rz(-1.8471142) q[1];
sx q[1];
rz(1.4599962) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0943992) q[0];
sx q[0];
rz(-1.5051821) q[0];
sx q[0];
rz(0.53334566) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47414549) q[2];
sx q[2];
rz(-0.80171004) q[2];
sx q[2];
rz(-2.2849883) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0467028) q[1];
sx q[1];
rz(-3.0071435) q[1];
sx q[1];
rz(1.0195169) q[1];
rz(3.0811062) q[3];
sx q[3];
rz(-2.2248259) q[3];
sx q[3];
rz(2.5791283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8755181) q[2];
sx q[2];
rz(-2.6660599) q[2];
sx q[2];
rz(-1.6391222) q[2];
rz(1.255704) q[3];
sx q[3];
rz(-1.4016822) q[3];
sx q[3];
rz(0.046253117) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2332377) q[0];
sx q[0];
rz(-0.73035705) q[0];
sx q[0];
rz(0.18049845) q[0];
rz(-1.3795229) q[1];
sx q[1];
rz(-2.5738398) q[1];
sx q[1];
rz(2.0464121) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4783673) q[0];
sx q[0];
rz(-2.992482) q[0];
sx q[0];
rz(-2.7115433) q[0];
x q[1];
rz(-2.344618) q[2];
sx q[2];
rz(-1.5124052) q[2];
sx q[2];
rz(-1.9820809) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.68183652) q[1];
sx q[1];
rz(-0.78826681) q[1];
sx q[1];
rz(-2.9440109) q[1];
x q[2];
rz(-0.88222701) q[3];
sx q[3];
rz(-1.467657) q[3];
sx q[3];
rz(2.8571199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3695662) q[2];
sx q[2];
rz(-0.84731805) q[2];
sx q[2];
rz(-2.0728716) q[2];
rz(-0.42701834) q[3];
sx q[3];
rz(-0.87978274) q[3];
sx q[3];
rz(-1.2515742) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1123947) q[0];
sx q[0];
rz(-0.91570941) q[0];
sx q[0];
rz(2.8503964) q[0];
rz(-1.4578106) q[1];
sx q[1];
rz(-1.0271415) q[1];
sx q[1];
rz(-0.88868946) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5989953) q[0];
sx q[0];
rz(-2.3610736) q[0];
sx q[0];
rz(0.88135834) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.027300731) q[2];
sx q[2];
rz(-0.75746398) q[2];
sx q[2];
rz(-0.94099076) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0093569) q[1];
sx q[1];
rz(-1.0601383) q[1];
sx q[1];
rz(0.64369323) q[1];
x q[2];
rz(0.8055775) q[3];
sx q[3];
rz(-0.42308261) q[3];
sx q[3];
rz(-2.7132636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7500744) q[2];
sx q[2];
rz(-1.3041648) q[2];
sx q[2];
rz(-2.1938426) q[2];
rz(-2.9446972) q[3];
sx q[3];
rz(-0.24053776) q[3];
sx q[3];
rz(-0.31015629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31203684) q[0];
sx q[0];
rz(-2.0289679) q[0];
sx q[0];
rz(2.6716676) q[0];
rz(-1.6795109) q[1];
sx q[1];
rz(-1.8406248) q[1];
sx q[1];
rz(-1.7376815) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3199661) q[0];
sx q[0];
rz(-1.2228773) q[0];
sx q[0];
rz(-2.1676284) q[0];
rz(-1.1833625) q[2];
sx q[2];
rz(-2.63113) q[2];
sx q[2];
rz(-1.1092157) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1822189) q[1];
sx q[1];
rz(-1.9314018) q[1];
sx q[1];
rz(-2.1350767) q[1];
rz(-0.49755712) q[3];
sx q[3];
rz(-0.64421875) q[3];
sx q[3];
rz(-1.8807172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.64138428) q[2];
sx q[2];
rz(-0.83157867) q[2];
sx q[2];
rz(-0.10759648) q[2];
rz(0.61947668) q[3];
sx q[3];
rz(-1.8457125) q[3];
sx q[3];
rz(1.7776325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3951931) q[0];
sx q[0];
rz(-1.0522319) q[0];
sx q[0];
rz(2.8885544) q[0];
rz(0.088317618) q[1];
sx q[1];
rz(-0.87366906) q[1];
sx q[1];
rz(-0.94892445) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.509453) q[0];
sx q[0];
rz(-1.4207736) q[0];
sx q[0];
rz(1.8709183) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1429092) q[2];
sx q[2];
rz(-2.6978931) q[2];
sx q[2];
rz(-2.2744409) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6151877) q[1];
sx q[1];
rz(-1.624409) q[1];
sx q[1];
rz(-2.0508906) q[1];
rz(2.1376417) q[3];
sx q[3];
rz(-1.1205871) q[3];
sx q[3];
rz(-0.70338827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39763149) q[2];
sx q[2];
rz(-2.2546015) q[2];
sx q[2];
rz(-2.8308947) q[2];
rz(-2.9001696) q[3];
sx q[3];
rz(-1.9390691) q[3];
sx q[3];
rz(-2.0588622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-0.13984118) q[0];
sx q[0];
rz(-2.7353291) q[0];
sx q[0];
rz(-1.9061506) q[0];
rz(-2.6605117) q[1];
sx q[1];
rz(-2.1886539) q[1];
sx q[1];
rz(-1.4280041) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9202703) q[0];
sx q[0];
rz(-2.3820665) q[0];
sx q[0];
rz(0.016245441) q[0];
rz(0.22869208) q[2];
sx q[2];
rz(-1.8426101) q[2];
sx q[2];
rz(0.52955176) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2216827) q[1];
sx q[1];
rz(-2.7836728) q[1];
sx q[1];
rz(2.329925) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72254027) q[3];
sx q[3];
rz(-1.3219537) q[3];
sx q[3];
rz(1.3470847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9674176) q[2];
sx q[2];
rz(-0.82505161) q[2];
sx q[2];
rz(0.027912557) q[2];
rz(2.2715955) q[3];
sx q[3];
rz(-2.3614707) q[3];
sx q[3];
rz(-1.1869259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1160527) q[0];
sx q[0];
rz(-1.5486131) q[0];
sx q[0];
rz(-1.0593587) q[0];
rz(1.4147883) q[1];
sx q[1];
rz(-1.6635514) q[1];
sx q[1];
rz(-0.47763225) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2842098) q[0];
sx q[0];
rz(-2.3071831) q[0];
sx q[0];
rz(-1.7718023) q[0];
x q[1];
rz(0.0019592802) q[2];
sx q[2];
rz(-1.9066487) q[2];
sx q[2];
rz(-2.4179469) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.76948386) q[1];
sx q[1];
rz(-1.4585351) q[1];
sx q[1];
rz(-3.0873507) q[1];
rz(0.60250256) q[3];
sx q[3];
rz(-0.91711603) q[3];
sx q[3];
rz(2.9127171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9719505) q[2];
sx q[2];
rz(-2.1596238) q[2];
sx q[2];
rz(2.7726445) q[2];
rz(1.2667027) q[3];
sx q[3];
rz(-2.4167175) q[3];
sx q[3];
rz(-0.96316159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13851588) q[0];
sx q[0];
rz(-1.930548) q[0];
sx q[0];
rz(1.7532274) q[0];
rz(-0.44454642) q[1];
sx q[1];
rz(-0.81041705) q[1];
sx q[1];
rz(3.0141426) q[1];
rz(-0.44345345) q[2];
sx q[2];
rz(-1.2724066) q[2];
sx q[2];
rz(2.3497697) q[2];
rz(-3.1034341) q[3];
sx q[3];
rz(-0.46783075) q[3];
sx q[3];
rz(-0.086188407) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
