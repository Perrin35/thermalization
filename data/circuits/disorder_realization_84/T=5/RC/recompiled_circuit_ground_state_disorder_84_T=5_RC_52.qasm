OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71275467) q[0];
sx q[0];
rz(-1.6102256) q[0];
sx q[0];
rz(1.1046326) q[0];
rz(1.334335) q[1];
sx q[1];
rz(-2.4185138) q[1];
sx q[1];
rz(-1.2637957) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0260213) q[0];
sx q[0];
rz(-1.8155671) q[0];
sx q[0];
rz(0.90306247) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99743263) q[2];
sx q[2];
rz(-1.9321529) q[2];
sx q[2];
rz(2.685355) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7967148) q[1];
sx q[1];
rz(-1.895322) q[1];
sx q[1];
rz(0.085436324) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24407152) q[3];
sx q[3];
rz(-1.9979949) q[3];
sx q[3];
rz(2.0352767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.15452142) q[2];
sx q[2];
rz(-2.0646586) q[2];
sx q[2];
rz(1.7729574) q[2];
rz(2.9030419) q[3];
sx q[3];
rz(-2.7261901) q[3];
sx q[3];
rz(-3.0273048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7411165) q[0];
sx q[0];
rz(-2.0023161) q[0];
sx q[0];
rz(1.9912632) q[0];
rz(-0.087609619) q[1];
sx q[1];
rz(-1.8116415) q[1];
sx q[1];
rz(1.5709343) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51975361) q[0];
sx q[0];
rz(-1.6770419) q[0];
sx q[0];
rz(0.48522093) q[0];
x q[1];
rz(1.5074128) q[2];
sx q[2];
rz(-1.3923613) q[2];
sx q[2];
rz(-3.1270535) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1521235) q[1];
sx q[1];
rz(-1.996576) q[1];
sx q[1];
rz(-2.2683737) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47111311) q[3];
sx q[3];
rz(-0.84976746) q[3];
sx q[3];
rz(1.6278933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.34438434) q[2];
sx q[2];
rz(-2.4281561) q[2];
sx q[2];
rz(-0.1864645) q[2];
rz(0.79948419) q[3];
sx q[3];
rz(-1.5954285) q[3];
sx q[3];
rz(2.1605087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28476533) q[0];
sx q[0];
rz(-1.3815877) q[0];
sx q[0];
rz(0.10936603) q[0];
rz(0.60375396) q[1];
sx q[1];
rz(-1.3394638) q[1];
sx q[1];
rz(2.107479) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4277735) q[0];
sx q[0];
rz(-1.7805011) q[0];
sx q[0];
rz(-0.13596491) q[0];
rz(1.6414406) q[2];
sx q[2];
rz(-1.6598668) q[2];
sx q[2];
rz(1.8235109) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.8960524) q[1];
sx q[1];
rz(-2.8421092) q[1];
sx q[1];
rz(-0.36548945) q[1];
rz(2.3699371) q[3];
sx q[3];
rz(-1.6724186) q[3];
sx q[3];
rz(1.1751428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5321396) q[2];
sx q[2];
rz(-0.9708465) q[2];
sx q[2];
rz(-2.5965221) q[2];
rz(-0.80101454) q[3];
sx q[3];
rz(-0.89371926) q[3];
sx q[3];
rz(-0.092122294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6545749) q[0];
sx q[0];
rz(-1.4734522) q[0];
sx q[0];
rz(2.6825478) q[0];
rz(-2.0252939) q[1];
sx q[1];
rz(-2.3479159) q[1];
sx q[1];
rz(-2.3150516) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5640663) q[0];
sx q[0];
rz(-0.56273976) q[0];
sx q[0];
rz(2.2731645) q[0];
x q[1];
rz(1.4722093) q[2];
sx q[2];
rz(-1.667983) q[2];
sx q[2];
rz(0.6310542) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.10545687) q[1];
sx q[1];
rz(-2.1107657) q[1];
sx q[1];
rz(2.9122874) q[1];
rz(-pi) q[2];
rz(2.8119254) q[3];
sx q[3];
rz(-0.35209823) q[3];
sx q[3];
rz(0.47131495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7894342) q[2];
sx q[2];
rz(-2.9106079) q[2];
sx q[2];
rz(-0.74971548) q[2];
rz(-2.0189144) q[3];
sx q[3];
rz(-1.9176982) q[3];
sx q[3];
rz(0.96021715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66626755) q[0];
sx q[0];
rz(-0.98639494) q[0];
sx q[0];
rz(3.0295897) q[0];
rz(-2.9169967) q[1];
sx q[1];
rz(-1.9738395) q[1];
sx q[1];
rz(0.78132838) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8664426) q[0];
sx q[0];
rz(-2.7366182) q[0];
sx q[0];
rz(1.5030131) q[0];
rz(-pi) q[1];
rz(-0.33022837) q[2];
sx q[2];
rz(-1.6598741) q[2];
sx q[2];
rz(1.1116127) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2558375) q[1];
sx q[1];
rz(-1.1449877) q[1];
sx q[1];
rz(0.7748697) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3329685) q[3];
sx q[3];
rz(-2.5356511) q[3];
sx q[3];
rz(1.1693418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3096699) q[2];
sx q[2];
rz(-2.2396542) q[2];
sx q[2];
rz(0.93264467) q[2];
rz(2.0617088) q[3];
sx q[3];
rz(-0.44411689) q[3];
sx q[3];
rz(3.1058969) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90750736) q[0];
sx q[0];
rz(-1.1312753) q[0];
sx q[0];
rz(1.750741) q[0];
rz(2.8098409) q[1];
sx q[1];
rz(-1.876588) q[1];
sx q[1];
rz(1.5501685) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0557779) q[0];
sx q[0];
rz(-1.5356488) q[0];
sx q[0];
rz(-2.2136392) q[0];
x q[1];
rz(2.6820002) q[2];
sx q[2];
rz(-1.6038648) q[2];
sx q[2];
rz(0.0065553105) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7227639) q[1];
sx q[1];
rz(-2.9420442) q[1];
sx q[1];
rz(0.90252374) q[1];
rz(-pi) q[2];
x q[2];
rz(2.329439) q[3];
sx q[3];
rz(-2.6366173) q[3];
sx q[3];
rz(2.2157089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8288237) q[2];
sx q[2];
rz(-2.2569816) q[2];
sx q[2];
rz(1.620232) q[2];
rz(2.17365) q[3];
sx q[3];
rz(-1.5827551) q[3];
sx q[3];
rz(-2.2255285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62436002) q[0];
sx q[0];
rz(-2.7150798) q[0];
sx q[0];
rz(-0.17159167) q[0];
rz(1.0393556) q[1];
sx q[1];
rz(-1.8828705) q[1];
sx q[1];
rz(-1.4412057) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87712578) q[0];
sx q[0];
rz(-1.4466337) q[0];
sx q[0];
rz(2.8178701) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51596291) q[2];
sx q[2];
rz(-1.4733286) q[2];
sx q[2];
rz(-0.24701842) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4644147) q[1];
sx q[1];
rz(-2.5485793) q[1];
sx q[1];
rz(0.89927425) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0073184) q[3];
sx q[3];
rz(-1.6905367) q[3];
sx q[3];
rz(-0.086825018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.27199304) q[2];
sx q[2];
rz(-1.4078434) q[2];
sx q[2];
rz(2.5970411) q[2];
rz(-2.6416687) q[3];
sx q[3];
rz(-2.1130424) q[3];
sx q[3];
rz(-2.669529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2917824) q[0];
sx q[0];
rz(-0.53684679) q[0];
sx q[0];
rz(-2.612402) q[0];
rz(-2.5657907) q[1];
sx q[1];
rz(-1.878783) q[1];
sx q[1];
rz(-2.0226488) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75615785) q[0];
sx q[0];
rz(-1.5917814) q[0];
sx q[0];
rz(0.051647112) q[0];
x q[1];
rz(0.37092692) q[2];
sx q[2];
rz(-1.7417522) q[2];
sx q[2];
rz(1.8425187) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0487489) q[1];
sx q[1];
rz(-2.1511937) q[1];
sx q[1];
rz(-1.012085) q[1];
x q[2];
rz(1.2203477) q[3];
sx q[3];
rz(-0.5548889) q[3];
sx q[3];
rz(-0.12400907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.90073663) q[2];
sx q[2];
rz(-2.6764328) q[2];
sx q[2];
rz(2.4294803) q[2];
rz(-1.5322878) q[3];
sx q[3];
rz(-2.1963547) q[3];
sx q[3];
rz(-1.7942662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3466472) q[0];
sx q[0];
rz(-2.0136588) q[0];
sx q[0];
rz(-1.0711063) q[0];
rz(1.935293) q[1];
sx q[1];
rz(-1.0164398) q[1];
sx q[1];
rz(1.4869022) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1893715) q[0];
sx q[0];
rz(-0.69714386) q[0];
sx q[0];
rz(-0.29336648) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63196147) q[2];
sx q[2];
rz(-2.2990531) q[2];
sx q[2];
rz(1.2974844) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2675954) q[1];
sx q[1];
rz(-1.2800299) q[1];
sx q[1];
rz(1.7520755) q[1];
rz(-pi) q[2];
rz(-2.1110299) q[3];
sx q[3];
rz(-2.6298454) q[3];
sx q[3];
rz(0.47357163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6335166) q[2];
sx q[2];
rz(-0.8997007) q[2];
sx q[2];
rz(0.077795204) q[2];
rz(-0.22732321) q[3];
sx q[3];
rz(-2.0096171) q[3];
sx q[3];
rz(1.0556861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8530497) q[0];
sx q[0];
rz(-2.396614) q[0];
sx q[0];
rz(-2.7689834) q[0];
rz(0.44652069) q[1];
sx q[1];
rz(-1.4563072) q[1];
sx q[1];
rz(-0.85652295) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3055182) q[0];
sx q[0];
rz(-1.5139607) q[0];
sx q[0];
rz(1.6147805) q[0];
rz(-2.2394453) q[2];
sx q[2];
rz(-2.9716316) q[2];
sx q[2];
rz(0.68357498) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0792993) q[1];
sx q[1];
rz(-0.41626272) q[1];
sx q[1];
rz(-1.8658616) q[1];
rz(-pi) q[2];
rz(-0.66121812) q[3];
sx q[3];
rz(-0.91877979) q[3];
sx q[3];
rz(2.2673502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7591758) q[2];
sx q[2];
rz(-2.0903812) q[2];
sx q[2];
rz(-2.5860533) q[2];
rz(-0.38604745) q[3];
sx q[3];
rz(-1.6470393) q[3];
sx q[3];
rz(2.4773795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1356708) q[0];
sx q[0];
rz(-1.4501403) q[0];
sx q[0];
rz(-1.8474664) q[0];
rz(0.80264965) q[1];
sx q[1];
rz(-1.4603271) q[1];
sx q[1];
rz(-2.7535798) q[1];
rz(1.6907202) q[2];
sx q[2];
rz(-1.9213866) q[2];
sx q[2];
rz(0.41768597) q[2];
rz(-0.84011806) q[3];
sx q[3];
rz(-1.5185322) q[3];
sx q[3];
rz(2.3867859) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
