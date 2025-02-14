OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.90091997) q[0];
sx q[0];
rz(-0.14578851) q[0];
sx q[0];
rz(1.2989651) q[0];
rz(2.9453912) q[1];
sx q[1];
rz(-1.8008404) q[1];
sx q[1];
rz(0.10032108) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.318905) q[0];
sx q[0];
rz(-1.5189369) q[0];
sx q[0];
rz(2.479631) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47122987) q[2];
sx q[2];
rz(-1.8719851) q[2];
sx q[2];
rz(-2.9848841) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70182204) q[1];
sx q[1];
rz(-1.7793097) q[1];
sx q[1];
rz(0.28065248) q[1];
x q[2];
rz(1.8842949) q[3];
sx q[3];
rz(-1.0035721) q[3];
sx q[3];
rz(2.6424266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4002865) q[2];
sx q[2];
rz(-0.15179673) q[2];
sx q[2];
rz(2.3347704) q[2];
rz(2.412879) q[3];
sx q[3];
rz(-2.3829134) q[3];
sx q[3];
rz(-1.2046643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62419409) q[0];
sx q[0];
rz(-0.26573467) q[0];
sx q[0];
rz(0.30701315) q[0];
rz(-1.864805) q[1];
sx q[1];
rz(-2.0025608) q[1];
sx q[1];
rz(2.0397287) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4808896) q[0];
sx q[0];
rz(-1.0202978) q[0];
sx q[0];
rz(-2.9376229) q[0];
x q[1];
rz(-2.8028433) q[2];
sx q[2];
rz(-0.18998665) q[2];
sx q[2];
rz(0.39443406) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4408177) q[1];
sx q[1];
rz(-1.7478746) q[1];
sx q[1];
rz(-0.13756474) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3086433) q[3];
sx q[3];
rz(-0.9303588) q[3];
sx q[3];
rz(0.73329496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.028194204) q[2];
sx q[2];
rz(-2.5598309) q[2];
sx q[2];
rz(-0.096435189) q[2];
rz(0.15245572) q[3];
sx q[3];
rz(-1.6343445) q[3];
sx q[3];
rz(2.4436387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0518799) q[0];
sx q[0];
rz(-1.1002325) q[0];
sx q[0];
rz(2.7413947) q[0];
rz(-1.6290889) q[1];
sx q[1];
rz(-2.9632603) q[1];
sx q[1];
rz(-2.8834744) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5494875) q[0];
sx q[0];
rz(-1.7096115) q[0];
sx q[0];
rz(-1.3347244) q[0];
rz(0.67419184) q[2];
sx q[2];
rz(-0.12756824) q[2];
sx q[2];
rz(1.3889165) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26545721) q[1];
sx q[1];
rz(-1.5578414) q[1];
sx q[1];
rz(1.575768) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3985004) q[3];
sx q[3];
rz(-1.8086284) q[3];
sx q[3];
rz(-2.3034277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1199946) q[2];
sx q[2];
rz(-1.1979411) q[2];
sx q[2];
rz(3.1094587) q[2];
rz(2.8739127) q[3];
sx q[3];
rz(-1.5175502) q[3];
sx q[3];
rz(1.5816429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5716008) q[0];
sx q[0];
rz(-1.6331693) q[0];
sx q[0];
rz(3.0573523) q[0];
rz(-3.1056504) q[1];
sx q[1];
rz(-0.033999559) q[1];
sx q[1];
rz(0.34119225) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0628898) q[0];
sx q[0];
rz(-0.69403946) q[0];
sx q[0];
rz(0.41578038) q[0];
x q[1];
rz(-2.7643599) q[2];
sx q[2];
rz(-2.5209941) q[2];
sx q[2];
rz(-2.008174) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.90662876) q[1];
sx q[1];
rz(-2.2412657) q[1];
sx q[1];
rz(2.5690298) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1197508) q[3];
sx q[3];
rz(-1.0195135) q[3];
sx q[3];
rz(-1.0663084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.42698947) q[2];
sx q[2];
rz(-2.0689071) q[2];
sx q[2];
rz(-1.6061456) q[2];
rz(2.8793907) q[3];
sx q[3];
rz(-1.5235135) q[3];
sx q[3];
rz(-2.1867627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3839805) q[0];
sx q[0];
rz(-0.42344991) q[0];
sx q[0];
rz(-1.0773995) q[0];
rz(-2.7491838) q[1];
sx q[1];
rz(-3.0632186) q[1];
sx q[1];
rz(1.0307301) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0466246) q[0];
sx q[0];
rz(-2.4343581) q[0];
sx q[0];
rz(-2.4695702) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0723128) q[2];
sx q[2];
rz(-1.4118115) q[2];
sx q[2];
rz(-1.0637525) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.64068078) q[1];
sx q[1];
rz(-1.7927153) q[1];
sx q[1];
rz(2.9725916) q[1];
rz(-2.0822273) q[3];
sx q[3];
rz(-1.2687195) q[3];
sx q[3];
rz(1.6405288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.82421676) q[2];
sx q[2];
rz(-0.65418303) q[2];
sx q[2];
rz(0.83924323) q[2];
rz(-0.9134891) q[3];
sx q[3];
rz(-1.313611) q[3];
sx q[3];
rz(0.14348468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6827253) q[0];
sx q[0];
rz(-2.904628) q[0];
sx q[0];
rz(1.6656026) q[0];
rz(-2.7492375) q[1];
sx q[1];
rz(-2.0456435) q[1];
sx q[1];
rz(-0.59142339) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05706035) q[0];
sx q[0];
rz(-2.1520237) q[0];
sx q[0];
rz(-2.6956062) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8663581) q[2];
sx q[2];
rz(-1.7448336) q[2];
sx q[2];
rz(-1.7013719) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5352262) q[1];
sx q[1];
rz(-2.6312345) q[1];
sx q[1];
rz(1.4354857) q[1];
rz(-pi) q[2];
rz(0.37844946) q[3];
sx q[3];
rz(-2.2317118) q[3];
sx q[3];
rz(1.4405516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.300294) q[2];
sx q[2];
rz(-2.4755307) q[2];
sx q[2];
rz(-2.5692614) q[2];
rz(2.9193997) q[3];
sx q[3];
rz(-2.7100345) q[3];
sx q[3];
rz(-0.64479327) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7543024) q[0];
sx q[0];
rz(-3.001725) q[0];
sx q[0];
rz(-2.733316) q[0];
rz(0.73221842) q[1];
sx q[1];
rz(-3.0156942) q[1];
sx q[1];
rz(-2.8439723) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2985122) q[0];
sx q[0];
rz(-1.1491927) q[0];
sx q[0];
rz(-1.9527736) q[0];
rz(2.0655259) q[2];
sx q[2];
rz(-1.552993) q[2];
sx q[2];
rz(-0.36843637) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.630328) q[1];
sx q[1];
rz(-1.2300092) q[1];
sx q[1];
rz(2.9326669) q[1];
rz(2.4378889) q[3];
sx q[3];
rz(-1.6856442) q[3];
sx q[3];
rz(-2.5435257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.97062651) q[2];
sx q[2];
rz(-1.2622702) q[2];
sx q[2];
rz(0.69620281) q[2];
rz(-2.2512186) q[3];
sx q[3];
rz(-1.9647157) q[3];
sx q[3];
rz(1.4674998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9625229) q[0];
sx q[0];
rz(-0.028554976) q[0];
sx q[0];
rz(0.21275511) q[0];
rz(0.46956024) q[1];
sx q[1];
rz(-2.1766365) q[1];
sx q[1];
rz(2.3874217) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6109723) q[0];
sx q[0];
rz(-1.7585131) q[0];
sx q[0];
rz(1.2769481) q[0];
rz(-pi) q[1];
rz(-0.32525678) q[2];
sx q[2];
rz(-1.7647247) q[2];
sx q[2];
rz(-3.0162899) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2060125) q[1];
sx q[1];
rz(-1.6846865) q[1];
sx q[1];
rz(2.4920032) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3776618) q[3];
sx q[3];
rz(-0.82909938) q[3];
sx q[3];
rz(1.0111077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0065877) q[2];
sx q[2];
rz(-2.2392515) q[2];
sx q[2];
rz(-0.78224409) q[2];
rz(-1.6953281) q[3];
sx q[3];
rz(-2.5952314) q[3];
sx q[3];
rz(2.7915891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049147216) q[0];
sx q[0];
rz(-0.47098422) q[0];
sx q[0];
rz(0.96442047) q[0];
rz(-1.2760705) q[1];
sx q[1];
rz(-1.4141021) q[1];
sx q[1];
rz(-1.5020348) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14012155) q[0];
sx q[0];
rz(-1.7709416) q[0];
sx q[0];
rz(1.8076623) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48351863) q[2];
sx q[2];
rz(-2.3502825) q[2];
sx q[2];
rz(0.4936337) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6148541) q[1];
sx q[1];
rz(-2.6236218) q[1];
sx q[1];
rz(1.6171842) q[1];
rz(-pi) q[2];
rz(1.1427059) q[3];
sx q[3];
rz(-0.16166281) q[3];
sx q[3];
rz(-0.77199329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2532578) q[2];
sx q[2];
rz(-1.2829245) q[2];
sx q[2];
rz(0.9453195) q[2];
rz(-2.3156598) q[3];
sx q[3];
rz(-1.5086987) q[3];
sx q[3];
rz(-2.6065684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.076040529) q[0];
sx q[0];
rz(-1.1689508) q[0];
sx q[0];
rz(-0.8031351) q[0];
rz(1.5777292) q[1];
sx q[1];
rz(-1.4811265) q[1];
sx q[1];
rz(-0.28958431) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3428925) q[0];
sx q[0];
rz(-0.10627593) q[0];
sx q[0];
rz(1.4551839) q[0];
x q[1];
rz(-1.2425735) q[2];
sx q[2];
rz(-1.8716836) q[2];
sx q[2];
rz(-2.2102578) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3037422) q[1];
sx q[1];
rz(-1.2134064) q[1];
sx q[1];
rz(2.9010495) q[1];
rz(-1.9042468) q[3];
sx q[3];
rz(-1.8419918) q[3];
sx q[3];
rz(1.9681794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3760066) q[2];
sx q[2];
rz(-3.0209318) q[2];
sx q[2];
rz(-2.181459) q[2];
rz(-0.58297408) q[3];
sx q[3];
rz(-2.4834902) q[3];
sx q[3];
rz(0.8647024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4509907) q[0];
sx q[0];
rz(-1.7091746) q[0];
sx q[0];
rz(-1.5102392) q[0];
rz(-3.1008537) q[1];
sx q[1];
rz(-0.67650411) q[1];
sx q[1];
rz(0.13112851) q[1];
rz(-2.1292674) q[2];
sx q[2];
rz(-0.53999117) q[2];
sx q[2];
rz(-2.7593031) q[2];
rz(1.4996281) q[3];
sx q[3];
rz(-1.6215848) q[3];
sx q[3];
rz(3.0231089) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
