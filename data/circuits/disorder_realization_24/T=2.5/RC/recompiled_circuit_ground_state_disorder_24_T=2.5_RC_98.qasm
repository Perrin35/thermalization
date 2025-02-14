OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3396575) q[0];
sx q[0];
rz(2.2428089) q[0];
sx q[0];
rz(8.1221683) q[0];
rz(-3.0175735) q[1];
sx q[1];
rz(-1.7350585) q[1];
sx q[1];
rz(-1.5136493) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4731015) q[0];
sx q[0];
rz(-1.2364179) q[0];
sx q[0];
rz(-2.8294493) q[0];
x q[1];
rz(2.3816125) q[2];
sx q[2];
rz(-0.70356762) q[2];
sx q[2];
rz(-2.4576996) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.096949654) q[1];
sx q[1];
rz(-1.5162807) q[1];
sx q[1];
rz(1.2300371) q[1];
rz(1.0694703) q[3];
sx q[3];
rz(-0.68911229) q[3];
sx q[3];
rz(-2.1769325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0579494) q[2];
sx q[2];
rz(-1.3336072) q[2];
sx q[2];
rz(-2.7570214) q[2];
rz(2.7390076) q[3];
sx q[3];
rz(-1.4339002) q[3];
sx q[3];
rz(-2.3334077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3835417) q[0];
sx q[0];
rz(-1.533968) q[0];
sx q[0];
rz(-0.59355271) q[0];
rz(1.3213762) q[1];
sx q[1];
rz(-1.079419) q[1];
sx q[1];
rz(0.039332565) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6932019) q[0];
sx q[0];
rz(-1.1993919) q[0];
sx q[0];
rz(-2.0601963) q[0];
rz(2.5051475) q[2];
sx q[2];
rz(-1.6628254) q[2];
sx q[2];
rz(-0.68816371) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.63168282) q[1];
sx q[1];
rz(-1.0541774) q[1];
sx q[1];
rz(1.5138018) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54941191) q[3];
sx q[3];
rz(-1.5447541) q[3];
sx q[3];
rz(0.56492912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3010657) q[2];
sx q[2];
rz(-1.7037921) q[2];
sx q[2];
rz(2.5145516) q[2];
rz(0.4380694) q[3];
sx q[3];
rz(-1.4984683) q[3];
sx q[3];
rz(1.1367249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42191926) q[0];
sx q[0];
rz(-1.0872343) q[0];
sx q[0];
rz(1.3713974) q[0];
rz(2.5321391) q[1];
sx q[1];
rz(-2.2513159) q[1];
sx q[1];
rz(1.9166463) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1430496) q[0];
sx q[0];
rz(-1.666905) q[0];
sx q[0];
rz(-0.19226464) q[0];
rz(-1.0880427) q[2];
sx q[2];
rz(-2.4039589) q[2];
sx q[2];
rz(2.9575153) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4546734) q[1];
sx q[1];
rz(-2.6489566) q[1];
sx q[1];
rz(-2.3083645) q[1];
x q[2];
rz(2.2663349) q[3];
sx q[3];
rz(-1.4039543) q[3];
sx q[3];
rz(0.60740208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3043392) q[2];
sx q[2];
rz(-0.86696583) q[2];
sx q[2];
rz(2.2875817) q[2];
rz(0.52418661) q[3];
sx q[3];
rz(-1.5093404) q[3];
sx q[3];
rz(-0.9700276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9147341) q[0];
sx q[0];
rz(-2.9840042) q[0];
sx q[0];
rz(2.3478813) q[0];
rz(0.0018250068) q[1];
sx q[1];
rz(-2.5853214) q[1];
sx q[1];
rz(-0.8108286) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2398259) q[0];
sx q[0];
rz(-1.7525273) q[0];
sx q[0];
rz(-2.316409) q[0];
rz(0.89703441) q[2];
sx q[2];
rz(-1.22222) q[2];
sx q[2];
rz(-0.86212117) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.27426611) q[1];
sx q[1];
rz(-1.9073448) q[1];
sx q[1];
rz(2.1907268) q[1];
rz(-pi) q[2];
rz(1.1417029) q[3];
sx q[3];
rz(-1.4388814) q[3];
sx q[3];
rz(-2.9219081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.68373716) q[2];
sx q[2];
rz(-2.808414) q[2];
sx q[2];
rz(1.6443171) q[2];
rz(1.3582683) q[3];
sx q[3];
rz(-2.0446916) q[3];
sx q[3];
rz(-1.8339405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7233647) q[0];
sx q[0];
rz(-1.6809373) q[0];
sx q[0];
rz(-0.25890589) q[0];
rz(-0.12340165) q[1];
sx q[1];
rz(-0.63107189) q[1];
sx q[1];
rz(-2.6554328) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21307316) q[0];
sx q[0];
rz(-0.80312906) q[0];
sx q[0];
rz(0.87246877) q[0];
rz(-pi) q[1];
rz(2.8672567) q[2];
sx q[2];
rz(-2.6736709) q[2];
sx q[2];
rz(0.58375025) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7546995) q[1];
sx q[1];
rz(-0.4122977) q[1];
sx q[1];
rz(1.0846433) q[1];
rz(2.5833273) q[3];
sx q[3];
rz(-2.1915806) q[3];
sx q[3];
rz(3.0594247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.90998021) q[2];
sx q[2];
rz(-2.6284802) q[2];
sx q[2];
rz(1.7463589) q[2];
rz(0.95476556) q[3];
sx q[3];
rz(-1.747811) q[3];
sx q[3];
rz(-2.098293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2024277) q[0];
sx q[0];
rz(-1.9068149) q[0];
sx q[0];
rz(2.3811316) q[0];
rz(2.3072534) q[1];
sx q[1];
rz(-2.0400679) q[1];
sx q[1];
rz(2.0371425) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9568312) q[0];
sx q[0];
rz(-2.5470877) q[0];
sx q[0];
rz(1.9678161) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.274674) q[2];
sx q[2];
rz(-2.6353177) q[2];
sx q[2];
rz(2.3828854) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.901746) q[1];
sx q[1];
rz(-2.2759328) q[1];
sx q[1];
rz(0.88259952) q[1];
rz(-pi) q[2];
rz(-2.831131) q[3];
sx q[3];
rz(-2.1230773) q[3];
sx q[3];
rz(2.3776835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1793648) q[2];
sx q[2];
rz(-0.70415512) q[2];
sx q[2];
rz(-0.041570138) q[2];
rz(2.9465594) q[3];
sx q[3];
rz(-0.2427559) q[3];
sx q[3];
rz(-2.354505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.70260173) q[0];
sx q[0];
rz(-2.2881303) q[0];
sx q[0];
rz(-2.3882197) q[0];
rz(-2.0808749) q[1];
sx q[1];
rz(-2.3371688) q[1];
sx q[1];
rz(-2.9341968) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8327861) q[0];
sx q[0];
rz(-0.96319234) q[0];
sx q[0];
rz(-0.73143801) q[0];
rz(-2.676187) q[2];
sx q[2];
rz(-1.8307181) q[2];
sx q[2];
rz(0.22451071) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8589391) q[1];
sx q[1];
rz(-0.6047073) q[1];
sx q[1];
rz(2.7207123) q[1];
rz(-pi) q[2];
rz(2.9459729) q[3];
sx q[3];
rz(-3.1017711) q[3];
sx q[3];
rz(0.34035836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.14219001) q[2];
sx q[2];
rz(-1.9839857) q[2];
sx q[2];
rz(-3.0863777) q[2];
rz(2.2986872) q[3];
sx q[3];
rz(-0.75777811) q[3];
sx q[3];
rz(-2.260476) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.371599) q[0];
sx q[0];
rz(-0.70796767) q[0];
sx q[0];
rz(1.9644894) q[0];
rz(-2.2616995) q[1];
sx q[1];
rz(-0.33029193) q[1];
sx q[1];
rz(-0.060861977) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56604474) q[0];
sx q[0];
rz(-1.4466929) q[0];
sx q[0];
rz(0.095416165) q[0];
x q[1];
rz(-1.251077) q[2];
sx q[2];
rz(-2.0315315) q[2];
sx q[2];
rz(1.2416897) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5338075) q[1];
sx q[1];
rz(-0.58064156) q[1];
sx q[1];
rz(0.21368475) q[1];
rz(-pi) q[2];
rz(-1.8375299) q[3];
sx q[3];
rz(-2.8681381) q[3];
sx q[3];
rz(-2.1463565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.047478288) q[2];
sx q[2];
rz(-0.55314174) q[2];
sx q[2];
rz(1.4609569) q[2];
rz(-0.85584062) q[3];
sx q[3];
rz(-2.1688921) q[3];
sx q[3];
rz(0.87876764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0395373) q[0];
sx q[0];
rz(-0.85298959) q[0];
sx q[0];
rz(-0.0041740388) q[0];
rz(1.8904842) q[1];
sx q[1];
rz(-0.1336385) q[1];
sx q[1];
rz(-2.4600696) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.174356) q[0];
sx q[0];
rz(-0.78034329) q[0];
sx q[0];
rz(-1.1514949) q[0];
rz(-pi) q[1];
rz(-2.5326742) q[2];
sx q[2];
rz(-1.0002197) q[2];
sx q[2];
rz(1.0341687) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.2118069) q[1];
sx q[1];
rz(-0.34084596) q[1];
sx q[1];
rz(-2.5548359) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.881991) q[3];
sx q[3];
rz(-2.3025945) q[3];
sx q[3];
rz(-1.9669718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.556813) q[2];
sx q[2];
rz(-2.6656373) q[2];
sx q[2];
rz(-1.3947831) q[2];
rz(0.41633385) q[3];
sx q[3];
rz(-1.2952015) q[3];
sx q[3];
rz(-2.7487315) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3861179) q[0];
sx q[0];
rz(-0.32805726) q[0];
sx q[0];
rz(-2.1395444) q[0];
rz(-2.877032) q[1];
sx q[1];
rz(-1.325565) q[1];
sx q[1];
rz(-1.4904259) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0430944) q[0];
sx q[0];
rz(-1.7014456) q[0];
sx q[0];
rz(2.5976973) q[0];
rz(-pi) q[1];
rz(0.80825808) q[2];
sx q[2];
rz(-1.4604476) q[2];
sx q[2];
rz(-2.7791952) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.893951) q[1];
sx q[1];
rz(-3.000285) q[1];
sx q[1];
rz(2.2013478) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6308925) q[3];
sx q[3];
rz(-1.4790241) q[3];
sx q[3];
rz(1.0453781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1313021) q[2];
sx q[2];
rz(-1.4069858) q[2];
sx q[2];
rz(-1.4170125) q[2];
rz(-0.33306444) q[3];
sx q[3];
rz(-0.76996961) q[3];
sx q[3];
rz(1.7634348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78111108) q[0];
sx q[0];
rz(-1.2662553) q[0];
sx q[0];
rz(1.8929831) q[0];
rz(-0.026451182) q[1];
sx q[1];
rz(-1.5298264) q[1];
sx q[1];
rz(-1.5419921) q[1];
rz(2.4290647) q[2];
sx q[2];
rz(-1.2488974) q[2];
sx q[2];
rz(1.2792994) q[2];
rz(-1.8881486) q[3];
sx q[3];
rz(-2.1614634) q[3];
sx q[3];
rz(0.95017016) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
