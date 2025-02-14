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
rz(-0.44218818) q[0];
sx q[0];
rz(-1.1531885) q[0];
sx q[0];
rz(-0.12864223) q[0];
rz(1.6755942) q[1];
sx q[1];
rz(-2.0425551) q[1];
sx q[1];
rz(-2.0597982) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8171514) q[0];
sx q[0];
rz(-1.6015581) q[0];
sx q[0];
rz(-3.0769303) q[0];
x q[1];
rz(2.2954582) q[2];
sx q[2];
rz(-0.48389176) q[2];
sx q[2];
rz(-2.0563375) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1342519) q[1];
sx q[1];
rz(-1.6414343) q[1];
sx q[1];
rz(1.411648) q[1];
rz(-pi) q[2];
rz(0.094688926) q[3];
sx q[3];
rz(-1.296954) q[3];
sx q[3];
rz(-1.7750334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.094540207) q[2];
sx q[2];
rz(-0.69539842) q[2];
sx q[2];
rz(-2.409234) q[2];
rz(0.098946027) q[3];
sx q[3];
rz(-1.0354592) q[3];
sx q[3];
rz(1.8100479) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.052213) q[0];
sx q[0];
rz(-0.78240028) q[0];
sx q[0];
rz(-0.041444929) q[0];
rz(0.6913569) q[1];
sx q[1];
rz(-1.8212916) q[1];
sx q[1];
rz(0.88821214) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9112638) q[0];
sx q[0];
rz(-0.8260051) q[0];
sx q[0];
rz(2.2780096) q[0];
rz(-pi) q[1];
rz(-2.0574548) q[2];
sx q[2];
rz(-2.1900822) q[2];
sx q[2];
rz(-0.82702434) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8687369) q[1];
sx q[1];
rz(-1.6140198) q[1];
sx q[1];
rz(-2.8985913) q[1];
rz(0.4732186) q[3];
sx q[3];
rz(-2.5788973) q[3];
sx q[3];
rz(1.109888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9899675) q[2];
sx q[2];
rz(-1.7513559) q[2];
sx q[2];
rz(1.3803231) q[2];
rz(0.049792854) q[3];
sx q[3];
rz(-2.4536665) q[3];
sx q[3];
rz(0.5602347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.77408537) q[0];
sx q[0];
rz(-0.15223509) q[0];
sx q[0];
rz(1.2138858) q[0];
rz(1.8269352) q[1];
sx q[1];
rz(-2.1036802) q[1];
sx q[1];
rz(-1.5705869) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6187046) q[0];
sx q[0];
rz(-1.9490593) q[0];
sx q[0];
rz(1.8357651) q[0];
rz(-pi) q[1];
rz(2.7309787) q[2];
sx q[2];
rz(-1.4272235) q[2];
sx q[2];
rz(-0.5856572) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5017058) q[1];
sx q[1];
rz(-0.16777953) q[1];
sx q[1];
rz(-1.3659992) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5475471) q[3];
sx q[3];
rz(-1.0540773) q[3];
sx q[3];
rz(-1.7917716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7870002) q[2];
sx q[2];
rz(-2.3921693) q[2];
sx q[2];
rz(-2.836239) q[2];
rz(-0.8461771) q[3];
sx q[3];
rz(-1.2140707) q[3];
sx q[3];
rz(2.2727216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8531891) q[0];
sx q[0];
rz(-2.1222293) q[0];
sx q[0];
rz(0.97002059) q[0];
rz(2.8963529) q[1];
sx q[1];
rz(-2.0742564) q[1];
sx q[1];
rz(-1.5018357) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1872699) q[0];
sx q[0];
rz(-1.4434955) q[0];
sx q[0];
rz(0.15592798) q[0];
rz(-1.7443329) q[2];
sx q[2];
rz(-1.8029986) q[2];
sx q[2];
rz(0.47877889) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0027155) q[1];
sx q[1];
rz(-1.7600696) q[1];
sx q[1];
rz(1.4156962) q[1];
rz(-pi) q[2];
rz(1.8957696) q[3];
sx q[3];
rz(-2.2712436) q[3];
sx q[3];
rz(2.9440299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5331427) q[2];
sx q[2];
rz(-1.5641944) q[2];
sx q[2];
rz(-1.2987785) q[2];
rz(0.53457824) q[3];
sx q[3];
rz(-2.5375073) q[3];
sx q[3];
rz(0.66876137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088148549) q[0];
sx q[0];
rz(-1.7867333) q[0];
sx q[0];
rz(-0.97766367) q[0];
rz(-2.4560302) q[1];
sx q[1];
rz(-0.79427636) q[1];
sx q[1];
rz(-0.96644863) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9622877) q[0];
sx q[0];
rz(-2.0317269) q[0];
sx q[0];
rz(2.3910752) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9990997) q[2];
sx q[2];
rz(-2.097528) q[2];
sx q[2];
rz(-1.525804) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.492474) q[1];
sx q[1];
rz(-0.8706514) q[1];
sx q[1];
rz(-2.9597069) q[1];
rz(-pi) q[2];
rz(2.6011916) q[3];
sx q[3];
rz(-2.5196919) q[3];
sx q[3];
rz(-2.2657713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5886249) q[2];
sx q[2];
rz(-1.7872461) q[2];
sx q[2];
rz(-2.8589613) q[2];
rz(-2.1770554) q[3];
sx q[3];
rz(-1.5805565) q[3];
sx q[3];
rz(-2.7717822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10837567) q[0];
sx q[0];
rz(-0.41154698) q[0];
sx q[0];
rz(2.0630398) q[0];
rz(-2.318553) q[1];
sx q[1];
rz(-2.0185202) q[1];
sx q[1];
rz(-2.3413234) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9554515) q[0];
sx q[0];
rz(-1.5793043) q[0];
sx q[0];
rz(-1.5999937) q[0];
rz(-0.3089463) q[2];
sx q[2];
rz(-1.2703203) q[2];
sx q[2];
rz(2.7555675) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.23689589) q[1];
sx q[1];
rz(-2.176732) q[1];
sx q[1];
rz(-0.69463457) q[1];
rz(-pi) q[2];
rz(-2.0144281) q[3];
sx q[3];
rz(-2.1882957) q[3];
sx q[3];
rz(1.8007743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0018953) q[2];
sx q[2];
rz(-1.8997833) q[2];
sx q[2];
rz(3.042172) q[2];
rz(-0.56441489) q[3];
sx q[3];
rz(-1.5494917) q[3];
sx q[3];
rz(0.91993371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1426706) q[0];
sx q[0];
rz(-0.26695928) q[0];
sx q[0];
rz(-3.1166742) q[0];
rz(1.092356) q[1];
sx q[1];
rz(-1.9905636) q[1];
sx q[1];
rz(-2.5999462) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.233802) q[0];
sx q[0];
rz(-2.4176717) q[0];
sx q[0];
rz(-2.1320599) q[0];
rz(1.0850026) q[2];
sx q[2];
rz(-2.5143904) q[2];
sx q[2];
rz(2.1202587) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.147824) q[1];
sx q[1];
rz(-1.2429825) q[1];
sx q[1];
rz(2.9891564) q[1];
rz(0.96407594) q[3];
sx q[3];
rz(-1.4243444) q[3];
sx q[3];
rz(-1.3769384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3148552) q[2];
sx q[2];
rz(-1.9172226) q[2];
sx q[2];
rz(-1.1401736) q[2];
rz(1.4402116) q[3];
sx q[3];
rz(-0.39339104) q[3];
sx q[3];
rz(-2.9668729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2575271) q[0];
sx q[0];
rz(-1.9767569) q[0];
sx q[0];
rz(0.87673941) q[0];
rz(-0.17403099) q[1];
sx q[1];
rz(-1.0935676) q[1];
sx q[1];
rz(1.3678975) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68326658) q[0];
sx q[0];
rz(-1.2948827) q[0];
sx q[0];
rz(1.0066731) q[0];
x q[1];
rz(0.22761818) q[2];
sx q[2];
rz(-1.926486) q[2];
sx q[2];
rz(1.5700036) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2217283) q[1];
sx q[1];
rz(-2.113453) q[1];
sx q[1];
rz(-0.43696398) q[1];
rz(-pi) q[2];
rz(-1.4488683) q[3];
sx q[3];
rz(-1.6462012) q[3];
sx q[3];
rz(0.46368515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4725264) q[2];
sx q[2];
rz(-1.3269227) q[2];
sx q[2];
rz(0.18829045) q[2];
rz(-1.7752198) q[3];
sx q[3];
rz(-1.7732737) q[3];
sx q[3];
rz(-1.9395456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-0.61757225) q[0];
sx q[0];
rz(-0.13372788) q[0];
sx q[0];
rz(-2.4753841) q[0];
rz(-1.1353525) q[1];
sx q[1];
rz(-2.1609781) q[1];
sx q[1];
rz(-1.9816144) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9177303) q[0];
sx q[0];
rz(-2.1356815) q[0];
sx q[0];
rz(-2.3333972) q[0];
rz(-pi) q[1];
rz(1.7889889) q[2];
sx q[2];
rz(-0.82880965) q[2];
sx q[2];
rz(0.11889501) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6659703) q[1];
sx q[1];
rz(-0.95487528) q[1];
sx q[1];
rz(3.1321976) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0951252) q[3];
sx q[3];
rz(-2.5216649) q[3];
sx q[3];
rz(1.7377095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.74405115) q[2];
sx q[2];
rz(-0.56637374) q[2];
sx q[2];
rz(2.9446824) q[2];
rz(0.87013733) q[3];
sx q[3];
rz(-1.0715485) q[3];
sx q[3];
rz(-1.2675233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63721913) q[0];
sx q[0];
rz(-2.0688031) q[0];
sx q[0];
rz(0.037242591) q[0];
rz(2.0540909) q[1];
sx q[1];
rz(-2.0934413) q[1];
sx q[1];
rz(0.16669272) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.113362) q[0];
sx q[0];
rz(-2.2489002) q[0];
sx q[0];
rz(1.2249169) q[0];
rz(1.8718821) q[2];
sx q[2];
rz(-1.1043806) q[2];
sx q[2];
rz(-0.25638858) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6065258) q[1];
sx q[1];
rz(-1.6589612) q[1];
sx q[1];
rz(2.8463581) q[1];
rz(-0.14504542) q[3];
sx q[3];
rz(-1.5022455) q[3];
sx q[3];
rz(-1.1762432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.4166261) q[2];
sx q[2];
rz(-2.5278957) q[2];
sx q[2];
rz(1.7566768) q[2];
rz(0.40063217) q[3];
sx q[3];
rz(-1.2847565) q[3];
sx q[3];
rz(2.1304776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.8818125) q[0];
sx q[0];
rz(-1.1008982) q[0];
sx q[0];
rz(-0.63006054) q[0];
rz(-0.13314816) q[1];
sx q[1];
rz(-1.1590191) q[1];
sx q[1];
rz(-1.0551183) q[1];
rz(-1.761663) q[2];
sx q[2];
rz(-1.6089572) q[2];
sx q[2];
rz(-1.937494) q[2];
rz(3.0146928) q[3];
sx q[3];
rz(-2.3681233) q[3];
sx q[3];
rz(-2.0518377) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
