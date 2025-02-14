OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.33534240722656) q[0];
sx q[0];
rz(3.5080533345514) q[0];
sx q[0];
rz(8.96564970015689) q[0];
rz(2.48998093605042) q[1];
sx q[1];
rz(4.8764570077234) q[1];
sx q[1];
rz(9.19303697942897) q[1];
cx q[1],q[0];
rz(0.153298616409302) q[0];
sx q[0];
rz(3.49203384120996) q[0];
sx q[0];
rz(11.1972086191098) q[0];
rz(1.66480207443237) q[2];
sx q[2];
rz(5.73478332360322) q[2];
sx q[2];
rz(8.03839085101291) q[2];
cx q[2],q[1];
rz(0.333179503679276) q[1];
sx q[1];
rz(1.19293788273866) q[1];
sx q[1];
rz(8.51986107825443) q[1];
rz(-3.10005521774292) q[3];
sx q[3];
rz(4.37429896195466) q[3];
sx q[3];
rz(11.0529561996381) q[3];
cx q[3],q[2];
rz(-0.0906524658203125) q[2];
sx q[2];
rz(3.67636892397935) q[2];
sx q[2];
rz(7.5621160030286) q[2];
rz(2.25179624557495) q[3];
sx q[3];
rz(4.66414049466188) q[3];
sx q[3];
rz(9.76340792178317) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.585625290870667) q[0];
sx q[0];
rz(5.37396851380403) q[0];
sx q[0];
rz(10.3505958080213) q[0];
rz(1.06491184234619) q[1];
sx q[1];
rz(4.71871856053407) q[1];
sx q[1];
rz(9.51025660186216) q[1];
cx q[1],q[0];
rz(-0.315622836351395) q[0];
sx q[0];
rz(5.1022026856714) q[0];
sx q[0];
rz(9.42833601044632) q[0];
rz(-0.477712899446487) q[2];
sx q[2];
rz(6.9517230113321) q[2];
sx q[2];
rz(8.60272619723483) q[2];
cx q[2],q[1];
rz(-1.8958455324173) q[1];
sx q[1];
rz(1.80080440838868) q[1];
sx q[1];
rz(8.88188318013355) q[1];
rz(-0.883933663368225) q[3];
sx q[3];
rz(5.08449641068513) q[3];
sx q[3];
rz(10.0372727274816) q[3];
cx q[3],q[2];
rz(1.81587898731232) q[2];
sx q[2];
rz(4.62058702309663) q[2];
sx q[2];
rz(12.9530346155088) q[2];
rz(1.92468416690826) q[3];
sx q[3];
rz(2.79721257288987) q[3];
sx q[3];
rz(9.64482845961257) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.33064079284668) q[0];
sx q[0];
rz(3.51295942266519) q[0];
sx q[0];
rz(9.92183644174739) q[0];
rz(2.09923839569092) q[1];
sx q[1];
rz(4.08916136820848) q[1];
sx q[1];
rz(9.13013327716991) q[1];
cx q[1],q[0];
rz(-0.932685852050781) q[0];
sx q[0];
rz(3.48104769189889) q[0];
sx q[0];
rz(10.787360405914) q[0];
rz(1.80148041248322) q[2];
sx q[2];
rz(2.4896224458986) q[2];
sx q[2];
rz(10.4204054236333) q[2];
cx q[2],q[1];
rz(4.49724149703979) q[1];
sx q[1];
rz(2.35017666419084) q[1];
sx q[1];
rz(8.79472813605472) q[1];
rz(-0.208468869328499) q[3];
sx q[3];
rz(3.75276819069917) q[3];
sx q[3];
rz(10.5465284347455) q[3];
cx q[3],q[2];
rz(-1.35433876514435) q[2];
sx q[2];
rz(5.42227354844148) q[2];
sx q[2];
rz(10.9865016698758) q[2];
rz(1.07623696327209) q[3];
sx q[3];
rz(4.44381669362123) q[3];
sx q[3];
rz(10.3537522911946) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.82660663127899) q[0];
sx q[0];
rz(4.7569159587198) q[0];
sx q[0];
rz(9.65409595369502) q[0];
rz(1.63531053066254) q[1];
sx q[1];
rz(4.82515934308107) q[1];
sx q[1];
rz(9.36270528509423) q[1];
cx q[1],q[0];
rz(2.69392251968384) q[0];
sx q[0];
rz(5.26644101937348) q[0];
sx q[0];
rz(10.4036947846334) q[0];
rz(-0.606212198734283) q[2];
sx q[2];
rz(5.81704154809053) q[2];
sx q[2];
rz(9.49852621405526) q[2];
cx q[2],q[1];
rz(2.55312991142273) q[1];
sx q[1];
rz(2.38479176362092) q[1];
sx q[1];
rz(8.26415977477237) q[1];
rz(-0.0932014435529709) q[3];
sx q[3];
rz(4.91249206860597) q[3];
sx q[3];
rz(10.1515598058622) q[3];
cx q[3],q[2];
rz(0.0929923579096794) q[2];
sx q[2];
rz(4.0576190670305) q[2];
sx q[2];
rz(10.7356645822446) q[2];
rz(0.0382895730435848) q[3];
sx q[3];
rz(4.93072024186189) q[3];
sx q[3];
rz(9.79939559697315) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.494860917329788) q[0];
sx q[0];
rz(2.42023876507814) q[0];
sx q[0];
rz(8.53178141116306) q[0];
rz(2.11217403411865) q[1];
sx q[1];
rz(2.41183772881562) q[1];
sx q[1];
rz(9.32889009862348) q[1];
cx q[1],q[0];
rz(-2.36781477928162) q[0];
sx q[0];
rz(3.826588960486) q[0];
sx q[0];
rz(10.1166401863019) q[0];
rz(4.6419506072998) q[2];
sx q[2];
rz(4.80109325249726) q[2];
sx q[2];
rz(9.08374903201267) q[2];
cx q[2],q[1];
rz(0.446605622768402) q[1];
sx q[1];
rz(4.16111996968324) q[1];
sx q[1];
rz(8.99456605910465) q[1];
rz(0.372627258300781) q[3];
sx q[3];
rz(5.19807139237458) q[3];
sx q[3];
rz(7.81181297301456) q[3];
cx q[3],q[2];
rz(0.154784947633743) q[2];
sx q[2];
rz(4.53165033658082) q[2];
sx q[2];
rz(8.99930099248096) q[2];
rz(-0.422438770532608) q[3];
sx q[3];
rz(4.2463550885492) q[3];
sx q[3];
rz(8.92160955666705) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.72858893871307) q[0];
sx q[0];
rz(2.5094030221277) q[0];
sx q[0];
rz(9.54194004683896) q[0];
rz(1.49401843547821) q[1];
sx q[1];
rz(4.04386767943437) q[1];
sx q[1];
rz(10.4959406614225) q[1];
cx q[1],q[0];
rz(3.32656311988831) q[0];
sx q[0];
rz(4.4302386363321) q[0];
sx q[0];
rz(8.35790405272647) q[0];
rz(1.12724959850311) q[2];
sx q[2];
rz(5.69505182107026) q[2];
sx q[2];
rz(9.33880751430198) q[2];
cx q[2],q[1];
rz(0.907296717166901) q[1];
sx q[1];
rz(4.9375298341089) q[1];
sx q[1];
rz(8.72721848487064) q[1];
rz(-1.25322306156158) q[3];
sx q[3];
rz(4.55671337445314) q[3];
sx q[3];
rz(8.91910687684222) q[3];
cx q[3],q[2];
rz(-2.56932067871094) q[2];
sx q[2];
rz(3.71650758584077) q[2];
sx q[2];
rz(12.5268869161527) q[2];
rz(-0.0764005407691002) q[3];
sx q[3];
rz(4.31140127976472) q[3];
sx q[3];
rz(12.1129424333493) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.70306783914566) q[0];
sx q[0];
rz(2.3301329930597) q[0];
sx q[0];
rz(10.3751208543698) q[0];
rz(1.15146613121033) q[1];
sx q[1];
rz(4.67154291470582) q[1];
sx q[1];
rz(10.7323304176252) q[1];
cx q[1],q[0];
rz(2.96145439147949) q[0];
sx q[0];
rz(2.74511200388009) q[0];
sx q[0];
rz(8.82759246825382) q[0];
rz(3.21443581581116) q[2];
sx q[2];
rz(1.19556394417817) q[2];
sx q[2];
rz(7.19240925311252) q[2];
cx q[2],q[1];
rz(3.55458354949951) q[1];
sx q[1];
rz(3.68281600077684) q[1];
sx q[1];
rz(6.70238873957797) q[1];
rz(1.43641412258148) q[3];
sx q[3];
rz(5.31409040291841) q[3];
sx q[3];
rz(9.84733015894099) q[3];
cx q[3],q[2];
rz(0.820397078990936) q[2];
sx q[2];
rz(3.82168457110459) q[2];
sx q[2];
rz(12.0127007722776) q[2];
rz(0.695948302745819) q[3];
sx q[3];
rz(4.61408594449098) q[3];
sx q[3];
rz(10.8799789905469) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.119479075074196) q[0];
sx q[0];
rz(4.83115604718263) q[0];
sx q[0];
rz(9.73168442248508) q[0];
rz(-1.27629971504211) q[1];
sx q[1];
rz(3.60797587235505) q[1];
sx q[1];
rz(11.3254184484403) q[1];
cx q[1],q[0];
rz(-1.06904447078705) q[0];
sx q[0];
rz(2.66761511762673) q[0];
sx q[0];
rz(9.77220696806117) q[0];
rz(0.972879230976105) q[2];
sx q[2];
rz(5.0737328847223) q[2];
sx q[2];
rz(8.12815520762607) q[2];
cx q[2],q[1];
rz(-0.136823743581772) q[1];
sx q[1];
rz(4.52421036561067) q[1];
sx q[1];
rz(7.53934524058505) q[1];
rz(0.960821032524109) q[3];
sx q[3];
rz(1.90533593495423) q[3];
sx q[3];
rz(11.3887028455655) q[3];
cx q[3],q[2];
rz(0.791573107242584) q[2];
sx q[2];
rz(2.3363523205095) q[2];
sx q[2];
rz(10.7150847673337) q[2];
rz(0.681141197681427) q[3];
sx q[3];
rz(4.04174295266206) q[3];
sx q[3];
rz(7.92297909259006) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.275555849075317) q[0];
sx q[0];
rz(5.24661269982392) q[0];
sx q[0];
rz(10.5984676837842) q[0];
rz(0.587000787258148) q[1];
sx q[1];
rz(2.36001101334626) q[1];
sx q[1];
rz(9.34506967513963) q[1];
cx q[1],q[0];
rz(0.226868867874146) q[0];
sx q[0];
rz(5.84392419655854) q[0];
sx q[0];
rz(12.0115859270017) q[0];
rz(0.112206839025021) q[2];
sx q[2];
rz(4.49229541619355) q[2];
sx q[2];
rz(12.8561069726865) q[2];
cx q[2],q[1];
rz(1.6708128452301) q[1];
sx q[1];
rz(5.30998125870759) q[1];
sx q[1];
rz(10.125480747215) q[1];
rz(-0.296717822551727) q[3];
sx q[3];
rz(3.928595276671) q[3];
sx q[3];
rz(10.2662244796674) q[3];
cx q[3],q[2];
rz(-0.540326356887817) q[2];
sx q[2];
rz(4.37022212346131) q[2];
sx q[2];
rz(10.7586953401487) q[2];
rz(1.59090828895569) q[3];
sx q[3];
rz(4.15167156060273) q[3];
sx q[3];
rz(9.61345914601489) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.482256680727005) q[0];
sx q[0];
rz(4.72005954583222) q[0];
sx q[0];
rz(7.57369194029971) q[0];
rz(0.0365296490490437) q[1];
sx q[1];
rz(4.30595842202241) q[1];
sx q[1];
rz(11.5220703840177) q[1];
cx q[1],q[0];
rz(-0.142577633261681) q[0];
sx q[0];
rz(1.86850884755189) q[0];
sx q[0];
rz(8.97828743457004) q[0];
rz(0.674514770507813) q[2];
sx q[2];
rz(2.72836286027963) q[2];
sx q[2];
rz(9.20880766808196) q[2];
cx q[2],q[1];
rz(-2.27303051948547) q[1];
sx q[1];
rz(2.29244628747041) q[1];
sx q[1];
rz(13.93123290538) q[1];
rz(0.114785276353359) q[3];
sx q[3];
rz(5.03944578965242) q[3];
sx q[3];
rz(9.59488322436019) q[3];
cx q[3],q[2];
rz(1.22699391841888) q[2];
sx q[2];
rz(4.17280867894227) q[2];
sx q[2];
rz(10.749556875221) q[2];
rz(0.756571829319) q[3];
sx q[3];
rz(5.394919308024) q[3];
sx q[3];
rz(10.2752619743268) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.0546773299574852) q[0];
sx q[0];
rz(5.43852558930451) q[0];
sx q[0];
rz(9.20050578414603) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-0.264515787363052) q[1];
sx q[1];
rz(5.0010336955362) q[1];
sx q[1];
rz(9.03411636351749) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-4.44431829452515) q[2];
sx q[2];
rz(2.78231674631173) q[2];
sx q[2];
rz(14.5773057699124) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(1.56394219398499) q[3];
sx q[3];
rz(5.3760885318094) q[3];
sx q[3];
rz(11.372741675369) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
