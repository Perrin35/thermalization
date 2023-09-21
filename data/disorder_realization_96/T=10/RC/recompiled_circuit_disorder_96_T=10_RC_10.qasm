OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.049790073) q[0];
sx q[0];
rz(3.0135305) q[0];
sx q[0];
rz(11.749) q[0];
rz(-2.1583537) q[1];
sx q[1];
rz(-2.6020738) q[1];
sx q[1];
rz(-1.9411545) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5759597) q[0];
sx q[0];
rz(-0.48086777) q[0];
sx q[0];
rz(0.54310449) q[0];
x q[1];
rz(-3.009216) q[2];
sx q[2];
rz(-2.1059603) q[2];
sx q[2];
rz(1.3995427) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9125036) q[1];
sx q[1];
rz(-2.2664245) q[1];
sx q[1];
rz(1.9181262) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8957542) q[3];
sx q[3];
rz(-2.339139) q[3];
sx q[3];
rz(1.6368395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0212705) q[2];
sx q[2];
rz(-2.5734084) q[2];
sx q[2];
rz(-1.583064) q[2];
rz(0.99672404) q[3];
sx q[3];
rz(-2.6895027) q[3];
sx q[3];
rz(-0.42580095) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5834171) q[0];
sx q[0];
rz(-2.3893864) q[0];
sx q[0];
rz(-3.0875207) q[0];
rz(-1.9460829) q[1];
sx q[1];
rz(-1.0369438) q[1];
sx q[1];
rz(0.53584677) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065952259) q[0];
sx q[0];
rz(-1.6951121) q[0];
sx q[0];
rz(-2.1527704) q[0];
rz(-pi) q[1];
rz(-0.70948647) q[2];
sx q[2];
rz(-1.3777133) q[2];
sx q[2];
rz(-0.62765861) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0332196) q[1];
sx q[1];
rz(-0.48711005) q[1];
sx q[1];
rz(2.5750722) q[1];
x q[2];
rz(-0.032580839) q[3];
sx q[3];
rz(-1.4279281) q[3];
sx q[3];
rz(2.8829765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1198931) q[2];
sx q[2];
rz(-2.0844441) q[2];
sx q[2];
rz(-2.1832441) q[2];
rz(-3.0751394) q[3];
sx q[3];
rz(-1.5575912) q[3];
sx q[3];
rz(2.691793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.7217343) q[0];
sx q[0];
rz(-1.2092713) q[0];
sx q[0];
rz(-0.15047519) q[0];
rz(2.6843605) q[1];
sx q[1];
rz(-0.91232863) q[1];
sx q[1];
rz(3.1157852) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31608554) q[0];
sx q[0];
rz(-1.8206017) q[0];
sx q[0];
rz(3.120963) q[0];
rz(3.0518968) q[2];
sx q[2];
rz(-1.8811474) q[2];
sx q[2];
rz(2.2657564) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5379996) q[1];
sx q[1];
rz(-0.94121658) q[1];
sx q[1];
rz(-2.0894719) q[1];
rz(1.953655) q[3];
sx q[3];
rz(-1.6930876) q[3];
sx q[3];
rz(-2.0910636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1228483) q[2];
sx q[2];
rz(-2.7763425) q[2];
sx q[2];
rz(-2.5562111) q[2];
rz(-0.18150005) q[3];
sx q[3];
rz(-1.8094962) q[3];
sx q[3];
rz(-1.5649149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90081763) q[0];
sx q[0];
rz(-2.5188991) q[0];
sx q[0];
rz(2.9649819) q[0];
rz(2.2606842) q[1];
sx q[1];
rz(-1.0662339) q[1];
sx q[1];
rz(2.6054629) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6648383) q[0];
sx q[0];
rz(-1.9739445) q[0];
sx q[0];
rz(-0.82157764) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23303194) q[2];
sx q[2];
rz(-2.8985902) q[2];
sx q[2];
rz(-2.255893) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.94169468) q[1];
sx q[1];
rz(-0.75921339) q[1];
sx q[1];
rz(-0.56337507) q[1];
rz(1.8648151) q[3];
sx q[3];
rz(-0.91563581) q[3];
sx q[3];
rz(-2.7340207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6716016) q[2];
sx q[2];
rz(-1.7207547) q[2];
sx q[2];
rz(-2.0969351) q[2];
rz(-2.4345543) q[3];
sx q[3];
rz(-0.98393327) q[3];
sx q[3];
rz(2.7846591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2351284) q[0];
sx q[0];
rz(-1.3236073) q[0];
sx q[0];
rz(-2.0902324) q[0];
rz(-1.6479187) q[1];
sx q[1];
rz(-0.58902478) q[1];
sx q[1];
rz(-0.043118127) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9263822) q[0];
sx q[0];
rz(-0.97391093) q[0];
sx q[0];
rz(-2.9125288) q[0];
rz(-2.0448858) q[2];
sx q[2];
rz(-2.5708963) q[2];
sx q[2];
rz(-1.2793465) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4700714) q[1];
sx q[1];
rz(-2.7579691) q[1];
sx q[1];
rz(2.3684711) q[1];
x q[2];
rz(2.352036) q[3];
sx q[3];
rz(-1.6766747) q[3];
sx q[3];
rz(-1.7136128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0126426) q[2];
sx q[2];
rz(-0.49762112) q[2];
sx q[2];
rz(-2.7894003) q[2];
rz(0.59018618) q[3];
sx q[3];
rz(-0.47368172) q[3];
sx q[3];
rz(-2.5804856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5181638) q[0];
sx q[0];
rz(-1.2151027) q[0];
sx q[0];
rz(-1.9859001) q[0];
rz(-2.39134) q[1];
sx q[1];
rz(-2.2022088) q[1];
sx q[1];
rz(2.0828784) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1832702) q[0];
sx q[0];
rz(-1.1407307) q[0];
sx q[0];
rz(1.3413315) q[0];
rz(-pi) q[1];
rz(-1.9361587) q[2];
sx q[2];
rz(-2.0632671) q[2];
sx q[2];
rz(-1.9838711) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8844879) q[1];
sx q[1];
rz(-0.37933644) q[1];
sx q[1];
rz(-0.92909716) q[1];
rz(-pi) q[2];
rz(1.3771463) q[3];
sx q[3];
rz(-0.87246694) q[3];
sx q[3];
rz(-2.8805672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6713312) q[2];
sx q[2];
rz(-1.417421) q[2];
sx q[2];
rz(-1.9227825) q[2];
rz(-1.9865215) q[3];
sx q[3];
rz(-0.23854908) q[3];
sx q[3];
rz(-1.7003805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041615151) q[0];
sx q[0];
rz(-1.3168553) q[0];
sx q[0];
rz(-1.51145) q[0];
rz(-1.7639683) q[1];
sx q[1];
rz(-2.8306077) q[1];
sx q[1];
rz(2.2999433) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.323303) q[0];
sx q[0];
rz(-2.0879732) q[0];
sx q[0];
rz(1.8316168) q[0];
x q[1];
rz(-2.878703) q[2];
sx q[2];
rz(-0.6436231) q[2];
sx q[2];
rz(2.7123244) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0591653) q[1];
sx q[1];
rz(-1.4945684) q[1];
sx q[1];
rz(3.0928897) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21861403) q[3];
sx q[3];
rz(-2.300005) q[3];
sx q[3];
rz(-1.5515755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1402309) q[2];
sx q[2];
rz(-1.3683687) q[2];
sx q[2];
rz(-3.0916396) q[2];
rz(-0.66155457) q[3];
sx q[3];
rz(-0.52246061) q[3];
sx q[3];
rz(0.18930999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89896232) q[0];
sx q[0];
rz(-2.1147418) q[0];
sx q[0];
rz(-1.7393973) q[0];
rz(0.095480355) q[1];
sx q[1];
rz(-1.1663368) q[1];
sx q[1];
rz(-0.41762525) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8861038) q[0];
sx q[0];
rz(-0.68414738) q[0];
sx q[0];
rz(-0.72649254) q[0];
x q[1];
rz(-0.11469658) q[2];
sx q[2];
rz(-1.9343997) q[2];
sx q[2];
rz(-0.39322688) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3216746) q[1];
sx q[1];
rz(-0.770861) q[1];
sx q[1];
rz(-0.8701156) q[1];
rz(0.91248625) q[3];
sx q[3];
rz(-0.45966002) q[3];
sx q[3];
rz(-1.0115136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3461356) q[2];
sx q[2];
rz(-1.6980349) q[2];
sx q[2];
rz(2.0987089) q[2];
rz(0.67388326) q[3];
sx q[3];
rz(-1.4833114) q[3];
sx q[3];
rz(-0.90014443) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3089356) q[0];
sx q[0];
rz(-1.7333663) q[0];
sx q[0];
rz(-0.62966627) q[0];
rz(0.57485238) q[1];
sx q[1];
rz(-1.31153) q[1];
sx q[1];
rz(-2.1946857) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6553584) q[0];
sx q[0];
rz(-1.2709193) q[0];
sx q[0];
rz(-2.8006058) q[0];
rz(-pi) q[1];
rz(0.38893716) q[2];
sx q[2];
rz(-0.10483327) q[2];
sx q[2];
rz(0.36738415) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2448333) q[1];
sx q[1];
rz(-1.4866801) q[1];
sx q[1];
rz(0.33478488) q[1];
x q[2];
rz(-1.1831207) q[3];
sx q[3];
rz(-2.6689853) q[3];
sx q[3];
rz(-1.2848867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.56069121) q[2];
sx q[2];
rz(-2.0118482) q[2];
sx q[2];
rz(1.2488731) q[2];
rz(-0.71436626) q[3];
sx q[3];
rz(-1.8656105) q[3];
sx q[3];
rz(-2.8760288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0666075) q[0];
sx q[0];
rz(-0.62244901) q[0];
sx q[0];
rz(2.4865436) q[0];
rz(-2.24522) q[1];
sx q[1];
rz(-0.90677774) q[1];
sx q[1];
rz(-0.64430976) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6039186) q[0];
sx q[0];
rz(-0.18225056) q[0];
sx q[0];
rz(0.3084348) q[0];
rz(-pi) q[1];
rz(-2.9115453) q[2];
sx q[2];
rz(-2.2736079) q[2];
sx q[2];
rz(-2.7237797) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.88313738) q[1];
sx q[1];
rz(-1.2346134) q[1];
sx q[1];
rz(1.2698445) q[1];
x q[2];
rz(-0.44390042) q[3];
sx q[3];
rz(-2.0022941) q[3];
sx q[3];
rz(2.3865226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7908988) q[2];
sx q[2];
rz(-1.472298) q[2];
sx q[2];
rz(-0.59990668) q[2];
rz(2.24263) q[3];
sx q[3];
rz(-0.18342429) q[3];
sx q[3];
rz(-1.775734) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29466378) q[0];
sx q[0];
rz(-1.8142721) q[0];
sx q[0];
rz(-0.66692764) q[0];
rz(-2.9121493) q[1];
sx q[1];
rz(-0.89090092) q[1];
sx q[1];
rz(0.13577239) q[1];
rz(-3.0145666) q[2];
sx q[2];
rz(-0.95022485) q[2];
sx q[2];
rz(-2.7765204) q[2];
rz(-2.7502144) q[3];
sx q[3];
rz(-0.93422514) q[3];
sx q[3];
rz(1.422062) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];