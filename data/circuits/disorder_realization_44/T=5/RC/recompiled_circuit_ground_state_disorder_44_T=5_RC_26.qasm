OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7527591) q[0];
sx q[0];
rz(1.692481) q[0];
sx q[0];
rz(11.115885) q[0];
rz(0.9736355) q[1];
sx q[1];
rz(-1.7042301) q[1];
sx q[1];
rz(-0.91926423) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96199233) q[0];
sx q[0];
rz(-1.5878979) q[0];
sx q[0];
rz(0.03058612) q[0];
rz(1.4798574) q[2];
sx q[2];
rz(-1.4509861) q[2];
sx q[2];
rz(2.7594942) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31499681) q[1];
sx q[1];
rz(-0.8657786) q[1];
sx q[1];
rz(-0.12048529) q[1];
rz(-pi) q[2];
rz(-1.1000865) q[3];
sx q[3];
rz(-1.9506427) q[3];
sx q[3];
rz(1.0468515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8093402) q[2];
sx q[2];
rz(-1.4490178) q[2];
sx q[2];
rz(-1.7285041) q[2];
rz(-2.938802) q[3];
sx q[3];
rz(-1.3821802) q[3];
sx q[3];
rz(-0.10281674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8973812) q[0];
sx q[0];
rz(-1.0845217) q[0];
sx q[0];
rz(-0.71075034) q[0];
rz(-2.3731025) q[1];
sx q[1];
rz(-1.0667421) q[1];
sx q[1];
rz(-2.1220727) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1508038) q[0];
sx q[0];
rz(-1.5312563) q[0];
sx q[0];
rz(-0.44724748) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4065845) q[2];
sx q[2];
rz(-2.7052042) q[2];
sx q[2];
rz(1.9917038) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.92495698) q[1];
sx q[1];
rz(-1.6018189) q[1];
sx q[1];
rz(-1.9418632) q[1];
x q[2];
rz(3.0583303) q[3];
sx q[3];
rz(-1.0181277) q[3];
sx q[3];
rz(-3.0577615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.76356137) q[2];
sx q[2];
rz(-1.2445933) q[2];
sx q[2];
rz(1.2223318) q[2];
rz(1.9289121) q[3];
sx q[3];
rz(-2.1247532) q[3];
sx q[3];
rz(2.4299664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0843622) q[0];
sx q[0];
rz(-2.1570692) q[0];
sx q[0];
rz(1.9236176) q[0];
rz(-0.51586622) q[1];
sx q[1];
rz(-2.5936544) q[1];
sx q[1];
rz(2.2191494) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9003214) q[0];
sx q[0];
rz(-1.5301955) q[0];
sx q[0];
rz(0.057828219) q[0];
rz(-pi) q[1];
rz(2.5848021) q[2];
sx q[2];
rz(-2.3344731) q[2];
sx q[2];
rz(-3.0592164) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.45097414) q[1];
sx q[1];
rz(-0.96267525) q[1];
sx q[1];
rz(0.7091503) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7397142) q[3];
sx q[3];
rz(-1.2460105) q[3];
sx q[3];
rz(-1.6095773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.018365232) q[2];
sx q[2];
rz(-1.0845228) q[2];
sx q[2];
rz(1.3398735) q[2];
rz(2.5943622) q[3];
sx q[3];
rz(-2.7180143) q[3];
sx q[3];
rz(2.4303998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0860586) q[0];
sx q[0];
rz(-0.14500293) q[0];
sx q[0];
rz(0.15922971) q[0];
rz(3.1314462) q[1];
sx q[1];
rz(-2.1187014) q[1];
sx q[1];
rz(2.8841282) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5437357) q[0];
sx q[0];
rz(-2.275995) q[0];
sx q[0];
rz(1.0969093) q[0];
x q[1];
rz(3.0765216) q[2];
sx q[2];
rz(-0.79539585) q[2];
sx q[2];
rz(-2.3415023) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.83492571) q[1];
sx q[1];
rz(-1.67027) q[1];
sx q[1];
rz(-0.62831459) q[1];
rz(0.61473989) q[3];
sx q[3];
rz(-2.2343544) q[3];
sx q[3];
rz(0.033294423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3782392) q[2];
sx q[2];
rz(-1.8469609) q[2];
sx q[2];
rz(-0.47284687) q[2];
rz(-2.4371448) q[3];
sx q[3];
rz(-1.364578) q[3];
sx q[3];
rz(2.4956467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6351629) q[0];
sx q[0];
rz(-1.9671054) q[0];
sx q[0];
rz(-3.0761062) q[0];
rz(-2.7323515) q[1];
sx q[1];
rz(-1.1301273) q[1];
sx q[1];
rz(-1.4261036) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62140853) q[0];
sx q[0];
rz(-1.6637633) q[0];
sx q[0];
rz(2.9220102) q[0];
x q[1];
rz(-1.6597366) q[2];
sx q[2];
rz(-1.1834025) q[2];
sx q[2];
rz(2.873444) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.1864422) q[1];
sx q[1];
rz(-1.3480061) q[1];
sx q[1];
rz(0.28526116) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0992005) q[3];
sx q[3];
rz(-0.33337731) q[3];
sx q[3];
rz(0.34587814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.19491974) q[2];
sx q[2];
rz(-1.0871004) q[2];
sx q[2];
rz(-1.8348414) q[2];
rz(-1.170018) q[3];
sx q[3];
rz(-2.7368059) q[3];
sx q[3];
rz(-3.034333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65866798) q[0];
sx q[0];
rz(-1.1558477) q[0];
sx q[0];
rz(-2.7401127) q[0];
rz(-0.67277706) q[1];
sx q[1];
rz(-0.79877001) q[1];
sx q[1];
rz(-0.85404095) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6333165) q[0];
sx q[0];
rz(-0.22249732) q[0];
sx q[0];
rz(-1.359324) q[0];
x q[1];
rz(1.7207022) q[2];
sx q[2];
rz(-1.7542363) q[2];
sx q[2];
rz(-2.4114087) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6352977) q[1];
sx q[1];
rz(-0.98165252) q[1];
sx q[1];
rz(1.6904669) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5307114) q[3];
sx q[3];
rz(-1.9457091) q[3];
sx q[3];
rz(1.553211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4868769) q[2];
sx q[2];
rz(-0.87149039) q[2];
sx q[2];
rz(-1.1775449) q[2];
rz(-2.5901637) q[3];
sx q[3];
rz(-1.5827554) q[3];
sx q[3];
rz(0.83561713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42171445) q[0];
sx q[0];
rz(-2.8604909) q[0];
sx q[0];
rz(1.1623435) q[0];
rz(-1.1516736) q[1];
sx q[1];
rz(-1.3538227) q[1];
sx q[1];
rz(-0.91845671) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3716482) q[0];
sx q[0];
rz(-1.838883) q[0];
sx q[0];
rz(1.2275342) q[0];
rz(-pi) q[1];
x q[1];
rz(0.8753885) q[2];
sx q[2];
rz(-1.3686485) q[2];
sx q[2];
rz(3.1234158) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9388401) q[1];
sx q[1];
rz(-1.563464) q[1];
sx q[1];
rz(-1.5640902) q[1];
x q[2];
rz(-1.2270532) q[3];
sx q[3];
rz(-1.4101646) q[3];
sx q[3];
rz(0.70987046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.034417001) q[2];
sx q[2];
rz(-1.7566661) q[2];
sx q[2];
rz(0.50298634) q[2];
rz(0.77477396) q[3];
sx q[3];
rz(-2.6796902) q[3];
sx q[3];
rz(-1.3431965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69304943) q[0];
sx q[0];
rz(-0.22668426) q[0];
sx q[0];
rz(1.6402798) q[0];
rz(-2.8003108) q[1];
sx q[1];
rz(-1.4309859) q[1];
sx q[1];
rz(1.1192082) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78748465) q[0];
sx q[0];
rz(-1.9409927) q[0];
sx q[0];
rz(-2.5404929) q[0];
rz(-pi) q[1];
rz(-2.0518752) q[2];
sx q[2];
rz(-1.8029034) q[2];
sx q[2];
rz(-0.26929917) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4367366) q[1];
sx q[1];
rz(-0.67187998) q[1];
sx q[1];
rz(0.19499548) q[1];
rz(0.31501997) q[3];
sx q[3];
rz(-0.4128939) q[3];
sx q[3];
rz(-1.1668432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1559653) q[2];
sx q[2];
rz(-0.95981821) q[2];
sx q[2];
rz(0.36925527) q[2];
rz(0.30820942) q[3];
sx q[3];
rz(-0.9674558) q[3];
sx q[3];
rz(1.5306028) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2823328) q[0];
sx q[0];
rz(-1.8349324) q[0];
sx q[0];
rz(-2.7985213) q[0];
rz(1.169091) q[1];
sx q[1];
rz(-1.3645376) q[1];
sx q[1];
rz(-1.4287359) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37821445) q[0];
sx q[0];
rz(-1.175048) q[0];
sx q[0];
rz(-1.2984492) q[0];
rz(0.83115432) q[2];
sx q[2];
rz(-2.2546283) q[2];
sx q[2];
rz(-1.9076951) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.27738849) q[1];
sx q[1];
rz(-2.7608747) q[1];
sx q[1];
rz(-1.4053132) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6713624) q[3];
sx q[3];
rz(-0.93358002) q[3];
sx q[3];
rz(1.5772082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2009361) q[2];
sx q[2];
rz(-2.9328465) q[2];
sx q[2];
rz(-1.7500056) q[2];
rz(-0.40677795) q[3];
sx q[3];
rz(-1.6986366) q[3];
sx q[3];
rz(-0.87882915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10130356) q[0];
sx q[0];
rz(-0.60281301) q[0];
sx q[0];
rz(-1.783675) q[0];
rz(1.2376002) q[1];
sx q[1];
rz(-2.1289181) q[1];
sx q[1];
rz(0.79992574) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7163776) q[0];
sx q[0];
rz(-2.3122462) q[0];
sx q[0];
rz(0.31405507) q[0];
rz(-pi) q[1];
rz(-1.8736035) q[2];
sx q[2];
rz(-1.0410415) q[2];
sx q[2];
rz(-2.985266) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7042735) q[1];
sx q[1];
rz(-1.4156439) q[1];
sx q[1];
rz(0.67351933) q[1];
x q[2];
rz(-1.1198197) q[3];
sx q[3];
rz(-0.52996892) q[3];
sx q[3];
rz(-1.3254904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.37821975) q[2];
sx q[2];
rz(-1.7174481) q[2];
sx q[2];
rz(0.073089449) q[2];
rz(2.8857005) q[3];
sx q[3];
rz(-2.3292694) q[3];
sx q[3];
rz(1.488744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8681317) q[0];
sx q[0];
rz(-1.68597) q[0];
sx q[0];
rz(-1.2706533) q[0];
rz(-1.801626) q[1];
sx q[1];
rz(-1.3506964) q[1];
sx q[1];
rz(0.66257308) q[1];
rz(-1.5727829) q[2];
sx q[2];
rz(-2.28021) q[2];
sx q[2];
rz(1.8328666) q[2];
rz(1.003391) q[3];
sx q[3];
rz(-2.1324674) q[3];
sx q[3];
rz(2.4218925) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
