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
rz(-2.2263865) q[0];
sx q[0];
rz(-2.6829166) q[0];
sx q[0];
rz(-0.31804481) q[0];
rz(7.6492352) q[1];
sx q[1];
rz(0.69861689) q[1];
sx q[1];
rz(3.0787992) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2554277) q[0];
sx q[0];
rz(-2.0768099) q[0];
sx q[0];
rz(-2.5939328) q[0];
rz(2.096296) q[2];
sx q[2];
rz(-2.3901148) q[2];
sx q[2];
rz(-0.95096248) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.54431984) q[1];
sx q[1];
rz(-1.2922799) q[1];
sx q[1];
rz(-1.0611978) q[1];
rz(1.0002238) q[3];
sx q[3];
rz(-2.4834342) q[3];
sx q[3];
rz(0.47129813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.41844765) q[2];
sx q[2];
rz(-1.3157088) q[2];
sx q[2];
rz(2.1212228) q[2];
rz(-2.4127035) q[3];
sx q[3];
rz(-0.43292361) q[3];
sx q[3];
rz(1.6093904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51434022) q[0];
sx q[0];
rz(-1.9753375) q[0];
sx q[0];
rz(3.1033206) q[0];
rz(-0.12380869) q[1];
sx q[1];
rz(-0.47272155) q[1];
sx q[1];
rz(1.5708057) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4306612) q[0];
sx q[0];
rz(-1.5607906) q[0];
sx q[0];
rz(-1.5512933) q[0];
x q[1];
rz(1.471631) q[2];
sx q[2];
rz(-2.3395774) q[2];
sx q[2];
rz(0.85814171) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5396186) q[1];
sx q[1];
rz(-1.569773) q[1];
sx q[1];
rz(3.1411693) q[1];
rz(1.4498715) q[3];
sx q[3];
rz(-2.0837415) q[3];
sx q[3];
rz(-2.4642022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.52246419) q[2];
sx q[2];
rz(-1.8121441) q[2];
sx q[2];
rz(0.5788571) q[2];
rz(2.0959496) q[3];
sx q[3];
rz(-2.3561616) q[3];
sx q[3];
rz(-2.3845909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66913644) q[0];
sx q[0];
rz(-2.0643015) q[0];
sx q[0];
rz(-0.45397595) q[0];
rz(2.9767735) q[1];
sx q[1];
rz(-2.3655128) q[1];
sx q[1];
rz(-1.4705315) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1810581) q[0];
sx q[0];
rz(-1.9583324) q[0];
sx q[0];
rz(2.9563532) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0963832) q[2];
sx q[2];
rz(-1.282935) q[2];
sx q[2];
rz(2.9381816) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.20023055) q[1];
sx q[1];
rz(-0.75298568) q[1];
sx q[1];
rz(-0.45238564) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64162125) q[3];
sx q[3];
rz(-2.372916) q[3];
sx q[3];
rz(-0.25594433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7620324) q[2];
sx q[2];
rz(-1.4780059) q[2];
sx q[2];
rz(1.3564159) q[2];
rz(-2.6421269) q[3];
sx q[3];
rz(-1.5288589) q[3];
sx q[3];
rz(1.0511901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.22719638) q[0];
sx q[0];
rz(-0.90407404) q[0];
sx q[0];
rz(1.0585349) q[0];
rz(-1.9805485) q[1];
sx q[1];
rz(-2.5807022) q[1];
sx q[1];
rz(-3.0243691) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49699052) q[0];
sx q[0];
rz(-1.0287026) q[0];
sx q[0];
rz(-1.5250456) q[0];
x q[1];
rz(1.8294705) q[2];
sx q[2];
rz(-0.66159464) q[2];
sx q[2];
rz(-1.0115397) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.30018729) q[1];
sx q[1];
rz(-0.8148199) q[1];
sx q[1];
rz(-2.7542265) q[1];
rz(2.5042512) q[3];
sx q[3];
rz(-1.9064603) q[3];
sx q[3];
rz(-0.46245271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3566572) q[2];
sx q[2];
rz(-1.6931809) q[2];
sx q[2];
rz(1.3680722) q[2];
rz(-1.28537) q[3];
sx q[3];
rz(-1.1077489) q[3];
sx q[3];
rz(2.4135597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042498978) q[0];
sx q[0];
rz(-0.92785257) q[0];
sx q[0];
rz(-1.4121144) q[0];
rz(-2.1389351) q[1];
sx q[1];
rz(-0.24191562) q[1];
sx q[1];
rz(1.5826506) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6928685) q[0];
sx q[0];
rz(-2.4131615) q[0];
sx q[0];
rz(2.4850575) q[0];
x q[1];
rz(-1.6539668) q[2];
sx q[2];
rz(-2.2874333) q[2];
sx q[2];
rz(0.73964529) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6465368) q[1];
sx q[1];
rz(-0.55083129) q[1];
sx q[1];
rz(-1.6158094) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2491426) q[3];
sx q[3];
rz(-1.6603004) q[3];
sx q[3];
rz(-0.22150515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5530508) q[2];
sx q[2];
rz(-1.7091227) q[2];
sx q[2];
rz(0.39349619) q[2];
rz(-0.95997512) q[3];
sx q[3];
rz(-2.4656651) q[3];
sx q[3];
rz(0.967832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8680962) q[0];
sx q[0];
rz(-0.12729004) q[0];
sx q[0];
rz(-2.9582276) q[0];
rz(2.6496437) q[1];
sx q[1];
rz(-2.349647) q[1];
sx q[1];
rz(-1.9221745) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8293835) q[0];
sx q[0];
rz(-1.6512733) q[0];
sx q[0];
rz(3.0521554) q[0];
rz(1.9398795) q[2];
sx q[2];
rz(-0.08412349) q[2];
sx q[2];
rz(2.8540128) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8811984) q[1];
sx q[1];
rz(-1.2165842) q[1];
sx q[1];
rz(-2.7071035) q[1];
rz(-pi) q[2];
rz(2.0424834) q[3];
sx q[3];
rz(-1.6460643) q[3];
sx q[3];
rz(1.4510297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.1313974) q[2];
sx q[2];
rz(-1.4364028) q[2];
sx q[2];
rz(0.08237002) q[2];
rz(-2.2149337) q[3];
sx q[3];
rz(-2.4121273) q[3];
sx q[3];
rz(-1.6389219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6019186) q[0];
sx q[0];
rz(-2.813756) q[0];
sx q[0];
rz(2.4251921) q[0];
rz(-0.40388233) q[1];
sx q[1];
rz(-1.5733746) q[1];
sx q[1];
rz(-1.4366038) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2879031) q[0];
sx q[0];
rz(-2.7349964) q[0];
sx q[0];
rz(-1.2000183) q[0];
x q[1];
rz(0.5282497) q[2];
sx q[2];
rz(-1.6852323) q[2];
sx q[2];
rz(-2.6334762) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3644412) q[1];
sx q[1];
rz(-1.4713411) q[1];
sx q[1];
rz(2.742139) q[1];
rz(-pi) q[2];
rz(2.3327504) q[3];
sx q[3];
rz(-2.1385962) q[3];
sx q[3];
rz(-1.1274991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7776103) q[2];
sx q[2];
rz(-2.1668375) q[2];
sx q[2];
rz(-0.065356143) q[2];
rz(-0.86723793) q[3];
sx q[3];
rz(-1.5012274) q[3];
sx q[3];
rz(-1.3245827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48677483) q[0];
sx q[0];
rz(-1.7205394) q[0];
sx q[0];
rz(2.4105893) q[0];
rz(0.77516088) q[1];
sx q[1];
rz(-0.33439264) q[1];
sx q[1];
rz(-0.41466546) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0867942) q[0];
sx q[0];
rz(-0.86547422) q[0];
sx q[0];
rz(-0.84673015) q[0];
rz(0.53040217) q[2];
sx q[2];
rz(-1.7549522) q[2];
sx q[2];
rz(-0.14094409) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29293167) q[1];
sx q[1];
rz(-0.39589105) q[1];
sx q[1];
rz(0.69940523) q[1];
rz(-pi) q[2];
rz(1.4067134) q[3];
sx q[3];
rz(-1.7030431) q[3];
sx q[3];
rz(-1.924072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.23184648) q[2];
sx q[2];
rz(-2.6555588) q[2];
sx q[2];
rz(3.0492142) q[2];
rz(-0.77629027) q[3];
sx q[3];
rz(-1.4624566) q[3];
sx q[3];
rz(-0.20364729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48731503) q[0];
sx q[0];
rz(-1.646811) q[0];
sx q[0];
rz(3.1238632) q[0];
rz(-2.0954633) q[1];
sx q[1];
rz(-1.4132063) q[1];
sx q[1];
rz(2.0808992) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8720406) q[0];
sx q[0];
rz(-1.5020796) q[0];
sx q[0];
rz(1.5567354) q[0];
x q[1];
rz(-1.9477884) q[2];
sx q[2];
rz(-0.62219884) q[2];
sx q[2];
rz(-2.5925328) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.13421397) q[1];
sx q[1];
rz(-1.7430647) q[1];
sx q[1];
rz(-2.3399347) q[1];
rz(-pi) q[2];
rz(-2.6178016) q[3];
sx q[3];
rz(-1.6219211) q[3];
sx q[3];
rz(1.2569951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1549418) q[2];
sx q[2];
rz(-0.81601802) q[2];
sx q[2];
rz(-0.6443392) q[2];
rz(-2.2996357) q[3];
sx q[3];
rz(-1.5717477) q[3];
sx q[3];
rz(-0.90698609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7613206) q[0];
sx q[0];
rz(-0.43571061) q[0];
sx q[0];
rz(0.45183387) q[0];
rz(2.1322346) q[1];
sx q[1];
rz(-1.511907) q[1];
sx q[1];
rz(1.570328) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46635755) q[0];
sx q[0];
rz(-1.5769099) q[0];
sx q[0];
rz(0.17019043) q[0];
rz(-pi) q[1];
rz(1.794593) q[2];
sx q[2];
rz(-2.6562956) q[2];
sx q[2];
rz(1.8388302) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45598509) q[1];
sx q[1];
rz(-2.7311014) q[1];
sx q[1];
rz(-1.2430771) q[1];
rz(-2.4191946) q[3];
sx q[3];
rz(-1.5649475) q[3];
sx q[3];
rz(1.4129782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62341225) q[2];
sx q[2];
rz(-2.0698915) q[2];
sx q[2];
rz(1.5709467) q[2];
rz(0.10562854) q[3];
sx q[3];
rz(-0.94625866) q[3];
sx q[3];
rz(2.6251729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8138206) q[0];
sx q[0];
rz(-2.1229424) q[0];
sx q[0];
rz(0.13737296) q[0];
rz(2.9227921) q[1];
sx q[1];
rz(-1.7411502) q[1];
sx q[1];
rz(2.2314744) q[1];
rz(1.1280288) q[2];
sx q[2];
rz(-1.3793066) q[2];
sx q[2];
rz(-2.1366091) q[2];
rz(1.2050592) q[3];
sx q[3];
rz(-2.1954721) q[3];
sx q[3];
rz(1.1385067) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
