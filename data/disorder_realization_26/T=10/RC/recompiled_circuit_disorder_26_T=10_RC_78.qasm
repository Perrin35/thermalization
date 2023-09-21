OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3044843) q[0];
sx q[0];
rz(-1.6882856) q[0];
sx q[0];
rz(2.8300571) q[0];
rz(-0.43752924) q[1];
sx q[1];
rz(-1.8234) q[1];
sx q[1];
rz(0.55895609) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1564724) q[0];
sx q[0];
rz(-0.44056842) q[0];
sx q[0];
rz(-1.23929) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5833601) q[2];
sx q[2];
rz(-2.1280834) q[2];
sx q[2];
rz(1.81665) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6701339) q[1];
sx q[1];
rz(-2.11073) q[1];
sx q[1];
rz(2.259841) q[1];
x q[2];
rz(2.7028014) q[3];
sx q[3];
rz(-1.3109129) q[3];
sx q[3];
rz(-2.4364542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3922334) q[2];
sx q[2];
rz(-1.8584676) q[2];
sx q[2];
rz(-2.5048845) q[2];
rz(-2.2926245) q[3];
sx q[3];
rz(-0.62148062) q[3];
sx q[3];
rz(-2.9076715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37671509) q[0];
sx q[0];
rz(-2.8945518) q[0];
sx q[0];
rz(0.15287457) q[0];
rz(2.3846467) q[1];
sx q[1];
rz(-1.5544954) q[1];
sx q[1];
rz(2.1551932) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61030412) q[0];
sx q[0];
rz(-3.0695519) q[0];
sx q[0];
rz(-0.59285592) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7631626) q[2];
sx q[2];
rz(-1.1726511) q[2];
sx q[2];
rz(-0.60000186) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0609329) q[1];
sx q[1];
rz(-2.6903209) q[1];
sx q[1];
rz(-2.5732451) q[1];
x q[2];
rz(3.0059079) q[3];
sx q[3];
rz(-1.8968582) q[3];
sx q[3];
rz(-2.5729834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5793005) q[2];
sx q[2];
rz(-1.9220756) q[2];
sx q[2];
rz(-2.3584649) q[2];
rz(0.018571818) q[3];
sx q[3];
rz(-1.5037856) q[3];
sx q[3];
rz(0.40772453) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0531533) q[0];
sx q[0];
rz(-2.8494371) q[0];
sx q[0];
rz(2.1799178) q[0];
rz(-0.36034521) q[1];
sx q[1];
rz(-1.1018437) q[1];
sx q[1];
rz(0.12869421) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012470915) q[0];
sx q[0];
rz(-3.0462628) q[0];
sx q[0];
rz(-1.0447797) q[0];
x q[1];
rz(-2.1305069) q[2];
sx q[2];
rz(-0.84583827) q[2];
sx q[2];
rz(-0.17979187) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8752746) q[1];
sx q[1];
rz(-0.53578636) q[1];
sx q[1];
rz(-2.5712625) q[1];
x q[2];
rz(0.80273654) q[3];
sx q[3];
rz(-1.667913) q[3];
sx q[3];
rz(0.10955284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0814357) q[2];
sx q[2];
rz(-1.1971985) q[2];
sx q[2];
rz(-1.241768) q[2];
rz(-2.5545819) q[3];
sx q[3];
rz(-0.95654064) q[3];
sx q[3];
rz(0.97755066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63242763) q[0];
sx q[0];
rz(-2.2594663) q[0];
sx q[0];
rz(1.084491) q[0];
rz(-1.658461) q[1];
sx q[1];
rz(-2.5741534) q[1];
sx q[1];
rz(-3.049057) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1922069) q[0];
sx q[0];
rz(-1.8398251) q[0];
sx q[0];
rz(-1.8640679) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7961568) q[2];
sx q[2];
rz(-2.0256809) q[2];
sx q[2];
rz(-2.5330184) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9738237) q[1];
sx q[1];
rz(-1.9134221) q[1];
sx q[1];
rz(-1.6325566) q[1];
rz(0.16102287) q[3];
sx q[3];
rz(-2.2194214) q[3];
sx q[3];
rz(-0.30740689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3080421) q[2];
sx q[2];
rz(-1.4414859) q[2];
sx q[2];
rz(-0.33205024) q[2];
rz(-2.0856693) q[3];
sx q[3];
rz(-0.27763405) q[3];
sx q[3];
rz(2.5312996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32245359) q[0];
sx q[0];
rz(-1.2912913) q[0];
sx q[0];
rz(2.9300368) q[0];
rz(-1.3062723) q[1];
sx q[1];
rz(-1.2439367) q[1];
sx q[1];
rz(-2.4938915) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0487329) q[0];
sx q[0];
rz(-1.6824241) q[0];
sx q[0];
rz(0.1603006) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34824246) q[2];
sx q[2];
rz(-1.7160545) q[2];
sx q[2];
rz(-2.9350304) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5375145) q[1];
sx q[1];
rz(-1.068183) q[1];
sx q[1];
rz(-0.22488774) q[1];
rz(-pi) q[2];
rz(-0.019142814) q[3];
sx q[3];
rz(-1.585841) q[3];
sx q[3];
rz(2.6263833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5806879) q[2];
sx q[2];
rz(-0.40955341) q[2];
sx q[2];
rz(2.4482751) q[2];
rz(0.66926113) q[3];
sx q[3];
rz(-1.7047434) q[3];
sx q[3];
rz(3.1366689) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82419056) q[0];
sx q[0];
rz(-0.22739246) q[0];
sx q[0];
rz(-1.90907) q[0];
rz(-1.0725853) q[1];
sx q[1];
rz(-1.0718081) q[1];
sx q[1];
rz(2.9673064) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.106819) q[0];
sx q[0];
rz(-0.79402059) q[0];
sx q[0];
rz(2.0110216) q[0];
x q[1];
rz(-0.66138791) q[2];
sx q[2];
rz(-0.97860133) q[2];
sx q[2];
rz(0.70450287) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7277158) q[1];
sx q[1];
rz(-1.1912279) q[1];
sx q[1];
rz(0.80645251) q[1];
rz(1.2586081) q[3];
sx q[3];
rz(-1.6815261) q[3];
sx q[3];
rz(2.1544416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.8217414) q[2];
sx q[2];
rz(-0.80703002) q[2];
sx q[2];
rz(2.9439587) q[2];
rz(-2.8526784) q[3];
sx q[3];
rz(-0.87696004) q[3];
sx q[3];
rz(-1.4060085) q[3];
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
rz(1.063868) q[0];
sx q[0];
rz(-2.6517695) q[0];
sx q[0];
rz(-0.20859627) q[0];
rz(-2.1754307) q[1];
sx q[1];
rz(-1.9816793) q[1];
sx q[1];
rz(1.6360412) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87058545) q[0];
sx q[0];
rz(-2.2858372) q[0];
sx q[0];
rz(-0.72512759) q[0];
x q[1];
rz(1.534329) q[2];
sx q[2];
rz(-2.0009396) q[2];
sx q[2];
rz(-1.9311116) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4475496) q[1];
sx q[1];
rz(-2.4195237) q[1];
sx q[1];
rz(0.85192792) q[1];
x q[2];
rz(-0.3397293) q[3];
sx q[3];
rz(-2.8390084) q[3];
sx q[3];
rz(1.8031977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2839526) q[2];
sx q[2];
rz(-0.74308926) q[2];
sx q[2];
rz(3.0440142) q[2];
rz(1.3939259) q[3];
sx q[3];
rz(-1.3503617) q[3];
sx q[3];
rz(-2.424749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7512648) q[0];
sx q[0];
rz(-1.3289691) q[0];
sx q[0];
rz(0.51399291) q[0];
rz(-3.0184074) q[1];
sx q[1];
rz(-2.8942278) q[1];
sx q[1];
rz(2.2095912) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65864627) q[0];
sx q[0];
rz(-2.6771149) q[0];
sx q[0];
rz(1.8770201) q[0];
rz(-pi) q[1];
rz(2.2064662) q[2];
sx q[2];
rz(-1.7889651) q[2];
sx q[2];
rz(0.93014923) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8646647) q[1];
sx q[1];
rz(-1.0229467) q[1];
sx q[1];
rz(0.95363708) q[1];
rz(2.2873015) q[3];
sx q[3];
rz(-1.7214509) q[3];
sx q[3];
rz(0.52779576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1121858) q[2];
sx q[2];
rz(-2.0307348) q[2];
sx q[2];
rz(-2.712148) q[2];
rz(-1.2094234) q[3];
sx q[3];
rz(-2.7691787) q[3];
sx q[3];
rz(-0.69303524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3747303) q[0];
sx q[0];
rz(-1.466789) q[0];
sx q[0];
rz(1.3388348) q[0];
rz(2.4354637) q[1];
sx q[1];
rz(-1.2495722) q[1];
sx q[1];
rz(-2.1910117) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75854077) q[0];
sx q[0];
rz(-1.0785111) q[0];
sx q[0];
rz(-2.5581215) q[0];
rz(-pi) q[1];
rz(-1.9865932) q[2];
sx q[2];
rz(-1.6396513) q[2];
sx q[2];
rz(1.5037675) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3705759) q[1];
sx q[1];
rz(-2.5713213) q[1];
sx q[1];
rz(-1.2194521) q[1];
x q[2];
rz(1.0334942) q[3];
sx q[3];
rz(-0.47035445) q[3];
sx q[3];
rz(-0.62894097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6616228) q[2];
sx q[2];
rz(-1.6500436) q[2];
sx q[2];
rz(-0.48428145) q[2];
rz(-2.2144923) q[3];
sx q[3];
rz(-1.8323332) q[3];
sx q[3];
rz(-0.62121975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1442239) q[0];
sx q[0];
rz(-0.08865083) q[0];
sx q[0];
rz(-0.22928672) q[0];
rz(0.43481049) q[1];
sx q[1];
rz(-1.2229342) q[1];
sx q[1];
rz(-0.71892175) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7867891) q[0];
sx q[0];
rz(-0.94593404) q[0];
sx q[0];
rz(0.11632365) q[0];
rz(-pi) q[1];
rz(-0.11833338) q[2];
sx q[2];
rz(-0.98630691) q[2];
sx q[2];
rz(-3.0362533) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.016991888) q[1];
sx q[1];
rz(-1.4233839) q[1];
sx q[1];
rz(-2.7684545) q[1];
x q[2];
rz(2.1925681) q[3];
sx q[3];
rz(-1.9927295) q[3];
sx q[3];
rz(-0.66872795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3832613) q[2];
sx q[2];
rz(-1.9394082) q[2];
sx q[2];
rz(0.74679217) q[2];
rz(2.2693999) q[3];
sx q[3];
rz(-2.3132497) q[3];
sx q[3];
rz(0.94223589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9185716) q[0];
sx q[0];
rz(-1.7445607) q[0];
sx q[0];
rz(1.8013409) q[0];
rz(2.7643798) q[1];
sx q[1];
rz(-1.4706392) q[1];
sx q[1];
rz(-2.6249862) q[1];
rz(1.8833075) q[2];
sx q[2];
rz(-1.0733114) q[2];
sx q[2];
rz(0.14807362) q[2];
rz(1.7198346) q[3];
sx q[3];
rz(-2.2825713) q[3];
sx q[3];
rz(1.5406516) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
