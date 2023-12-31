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
rz(-0.31153554) q[0];
rz(2.7040634) q[1];
sx q[1];
rz(-1.3181926) q[1];
sx q[1];
rz(-0.55895609) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9851202) q[0];
sx q[0];
rz(-2.7010242) q[0];
sx q[0];
rz(-1.23929) q[0];
x q[1];
rz(-1.5582325) q[2];
sx q[2];
rz(-2.1280834) q[2];
sx q[2];
rz(-1.81665) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.47145876) q[1];
sx q[1];
rz(-1.0308627) q[1];
sx q[1];
rz(-2.259841) q[1];
rz(0.43879126) q[3];
sx q[3];
rz(-1.3109129) q[3];
sx q[3];
rz(-0.70513844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7493593) q[2];
sx q[2];
rz(-1.2831251) q[2];
sx q[2];
rz(-2.5048845) q[2];
rz(-0.84896815) q[3];
sx q[3];
rz(-2.520112) q[3];
sx q[3];
rz(-2.9076715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7648776) q[0];
sx q[0];
rz(-0.24704084) q[0];
sx q[0];
rz(-0.15287457) q[0];
rz(0.75694594) q[1];
sx q[1];
rz(-1.5544954) q[1];
sx q[1];
rz(-2.1551932) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5312885) q[0];
sx q[0];
rz(-0.072040759) q[0];
sx q[0];
rz(-2.5487367) q[0];
x q[1];
rz(0.54665357) q[2];
sx q[2];
rz(-2.2998027) q[2];
sx q[2];
rz(0.58568776) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0806597) q[1];
sx q[1];
rz(-2.6903209) q[1];
sx q[1];
rz(-0.56834759) q[1];
x q[2];
rz(1.1902477) q[3];
sx q[3];
rz(-2.7893587) q[3];
sx q[3];
rz(-0.16570839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5622921) q[2];
sx q[2];
rz(-1.9220756) q[2];
sx q[2];
rz(2.3584649) q[2];
rz(-0.018571818) q[3];
sx q[3];
rz(-1.5037856) q[3];
sx q[3];
rz(-0.40772453) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0531533) q[0];
sx q[0];
rz(-0.29215559) q[0];
sx q[0];
rz(0.96167481) q[0];
rz(-0.36034521) q[1];
sx q[1];
rz(-2.0397489) q[1];
sx q[1];
rz(-0.12869421) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0823682) q[0];
sx q[0];
rz(-1.5229862) q[0];
sx q[0];
rz(-1.6533018) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0110858) q[2];
sx q[2];
rz(-2.2957544) q[2];
sx q[2];
rz(-2.9618008) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8752746) q[1];
sx q[1];
rz(-0.53578636) q[1];
sx q[1];
rz(-2.5712625) q[1];
x q[2];
rz(-1.7101173) q[3];
sx q[3];
rz(-0.77292597) q[3];
sx q[3];
rz(-1.7804002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0814357) q[2];
sx q[2];
rz(-1.9443941) q[2];
sx q[2];
rz(1.241768) q[2];
rz(-0.5870108) q[3];
sx q[3];
rz(-2.185052) q[3];
sx q[3];
rz(0.97755066) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63242763) q[0];
sx q[0];
rz(-0.88212633) q[0];
sx q[0];
rz(-2.0571016) q[0];
rz(-1.658461) q[1];
sx q[1];
rz(-0.56743923) q[1];
sx q[1];
rz(-0.09253563) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0412484) q[0];
sx q[0];
rz(-2.7462602) q[0];
sx q[0];
rz(-0.80907099) q[0];
rz(-pi) q[1];
rz(-2.0501577) q[2];
sx q[2];
rz(-1.2617246) q[2];
sx q[2];
rz(-1.1190337) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.167769) q[1];
sx q[1];
rz(-1.9134221) q[1];
sx q[1];
rz(-1.6325566) q[1];
x q[2];
rz(2.9805698) q[3];
sx q[3];
rz(-0.92217126) q[3];
sx q[3];
rz(2.8341858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83355054) q[2];
sx q[2];
rz(-1.7001067) q[2];
sx q[2];
rz(-0.33205024) q[2];
rz(1.0559233) q[3];
sx q[3];
rz(-0.27763405) q[3];
sx q[3];
rz(2.5312996) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-1.897656) q[1];
sx q[1];
rz(2.4938915) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0603795) q[0];
sx q[0];
rz(-0.19506422) q[0];
sx q[0];
rz(0.61225201) q[0];
x q[1];
rz(-2.7933502) q[2];
sx q[2];
rz(-1.7160545) q[2];
sx q[2];
rz(-0.20656221) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.065121) q[1];
sx q[1];
rz(-1.7674801) q[1];
sx q[1];
rz(-2.0842488) q[1];
rz(2.475431) q[3];
sx q[3];
rz(-3.1172459) q[3];
sx q[3];
rz(1.7216046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5806879) q[2];
sx q[2];
rz(-2.7320392) q[2];
sx q[2];
rz(-0.69331759) q[2];
rz(-2.4723315) q[3];
sx q[3];
rz(-1.7047434) q[3];
sx q[3];
rz(-0.0049237331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3174021) q[0];
sx q[0];
rz(-0.22739246) q[0];
sx q[0];
rz(-1.2325226) q[0];
rz(-1.0725853) q[1];
sx q[1];
rz(-1.0718081) q[1];
sx q[1];
rz(2.9673064) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3586853) q[0];
sx q[0];
rz(-1.2619962) q[0];
sx q[0];
rz(-2.3147644) q[0];
rz(-pi) q[1];
rz(-2.2767378) q[2];
sx q[2];
rz(-2.1055429) q[2];
sx q[2];
rz(-0.45644444) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4982521) q[1];
sx q[1];
rz(-2.2689515) q[1];
sx q[1];
rz(0.50486418) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2234736) q[3];
sx q[3];
rz(-2.8109549) q[3];
sx q[3];
rz(-0.91352458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.8217414) q[2];
sx q[2];
rz(-0.80703002) q[2];
sx q[2];
rz(-0.19763395) q[2];
rz(-0.28891426) q[3];
sx q[3];
rz(-0.87696004) q[3];
sx q[3];
rz(-1.7355841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.063868) q[0];
sx q[0];
rz(-2.6517695) q[0];
sx q[0];
rz(-2.9329964) q[0];
rz(0.96616191) q[1];
sx q[1];
rz(-1.1599133) q[1];
sx q[1];
rz(-1.6360412) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9150328) q[0];
sx q[0];
rz(-1.0462927) q[0];
sx q[0];
rz(-2.4302308) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43039544) q[2];
sx q[2];
rz(-1.5376523) q[2];
sx q[2];
rz(0.37552777) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.70430763) q[1];
sx q[1];
rz(-1.1204801) q[1];
sx q[1];
rz(-0.98547658) q[1];
rz(-pi) q[2];
rz(-0.286245) q[3];
sx q[3];
rz(-1.6702594) q[3];
sx q[3];
rz(-0.092982987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2839526) q[2];
sx q[2];
rz(-2.3985034) q[2];
sx q[2];
rz(0.097578438) q[2];
rz(-1.7476667) q[3];
sx q[3];
rz(-1.3503617) q[3];
sx q[3];
rz(0.71684366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7512648) q[0];
sx q[0];
rz(-1.8126235) q[0];
sx q[0];
rz(-2.6275997) q[0];
rz(0.12318525) q[1];
sx q[1];
rz(-2.8942278) q[1];
sx q[1];
rz(-0.93200144) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4829464) q[0];
sx q[0];
rz(-0.46447771) q[0];
sx q[0];
rz(-1.8770201) q[0];
rz(-pi) q[1];
rz(2.8727505) q[2];
sx q[2];
rz(-0.95250597) q[2];
sx q[2];
rz(0.48228574) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8646647) q[1];
sx q[1];
rz(-1.0229467) q[1];
sx q[1];
rz(-0.95363708) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85429116) q[3];
sx q[3];
rz(-1.7214509) q[3];
sx q[3];
rz(2.6137969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1121858) q[2];
sx q[2];
rz(-2.0307348) q[2];
sx q[2];
rz(-0.4294447) q[2];
rz(-1.9321692) q[3];
sx q[3];
rz(-0.37241396) q[3];
sx q[3];
rz(-0.69303524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76686239) q[0];
sx q[0];
rz(-1.466789) q[0];
sx q[0];
rz(-1.3388348) q[0];
rz(-2.4354637) q[1];
sx q[1];
rz(-1.2495722) q[1];
sx q[1];
rz(2.1910117) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3830519) q[0];
sx q[0];
rz(-1.0785111) q[0];
sx q[0];
rz(-2.5581215) q[0];
rz(1.1549994) q[2];
sx q[2];
rz(-1.6396513) q[2];
sx q[2];
rz(-1.6378251) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7710167) q[1];
sx q[1];
rz(-0.5702714) q[1];
sx q[1];
rz(1.9221406) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1589963) q[3];
sx q[3];
rz(-1.3367062) q[3];
sx q[3];
rz(1.4300508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6616228) q[2];
sx q[2];
rz(-1.491549) q[2];
sx q[2];
rz(0.48428145) q[2];
rz(-2.2144923) q[3];
sx q[3];
rz(-1.3092594) q[3];
sx q[3];
rz(0.62121975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.1442239) q[0];
sx q[0];
rz(-3.0529418) q[0];
sx q[0];
rz(2.9123059) q[0];
rz(-0.43481049) q[1];
sx q[1];
rz(-1.2229342) q[1];
sx q[1];
rz(0.71892175) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9839448) q[0];
sx q[0];
rz(-0.63417182) q[0];
sx q[0];
rz(1.7303403) q[0];
rz(-pi) q[1];
rz(-0.98307307) q[2];
sx q[2];
rz(-1.6694153) q[2];
sx q[2];
rz(1.3999511) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6112411) q[1];
sx q[1];
rz(-1.9396922) q[1];
sx q[1];
rz(-1.4126652) q[1];
x q[2];
rz(-0.9142011) q[3];
sx q[3];
rz(-0.73540348) q[3];
sx q[3];
rz(2.7587193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3832613) q[2];
sx q[2];
rz(-1.2021844) q[2];
sx q[2];
rz(2.3948005) q[2];
rz(0.87219277) q[3];
sx q[3];
rz(-0.82834297) q[3];
sx q[3];
rz(0.94223589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9185716) q[0];
sx q[0];
rz(-1.3970319) q[0];
sx q[0];
rz(-1.3402517) q[0];
rz(-2.7643798) q[1];
sx q[1];
rz(-1.6709534) q[1];
sx q[1];
rz(0.51660641) q[1];
rz(1.8833075) q[2];
sx q[2];
rz(-1.0733114) q[2];
sx q[2];
rz(0.14807362) q[2];
rz(1.4217581) q[3];
sx q[3];
rz(-0.85902135) q[3];
sx q[3];
rz(-1.600941) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
