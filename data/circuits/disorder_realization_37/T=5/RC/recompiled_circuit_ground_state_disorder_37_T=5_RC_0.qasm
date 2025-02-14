OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6310298) q[0];
sx q[0];
rz(-2.6369542) q[0];
sx q[0];
rz(-2.7161427) q[0];
rz(0.57234859) q[1];
sx q[1];
rz(-2.0444137) q[1];
sx q[1];
rz(-1.9414577) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10353032) q[0];
sx q[0];
rz(-1.5804884) q[0];
sx q[0];
rz(-0.25083019) q[0];
x q[1];
rz(-1.9169943) q[2];
sx q[2];
rz(-1.8270512) q[2];
sx q[2];
rz(3.0145149) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3878183) q[1];
sx q[1];
rz(-1.1987616) q[1];
sx q[1];
rz(-1.340999) q[1];
rz(-pi) q[2];
rz(0.17961924) q[3];
sx q[3];
rz(-2.4625375) q[3];
sx q[3];
rz(1.8906901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9878865) q[2];
sx q[2];
rz(-1.989863) q[2];
sx q[2];
rz(-2.494334) q[2];
rz(-0.17793947) q[3];
sx q[3];
rz(-1.2578332) q[3];
sx q[3];
rz(1.9378763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-1.7385638) q[0];
sx q[0];
rz(-0.084349923) q[0];
sx q[0];
rz(2.6179598) q[0];
rz(-1.9117484) q[1];
sx q[1];
rz(-2.3975027) q[1];
sx q[1];
rz(-2.8935166) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6726674) q[0];
sx q[0];
rz(-1.6771731) q[0];
sx q[0];
rz(-2.399978) q[0];
rz(-pi) q[1];
rz(-1.1311943) q[2];
sx q[2];
rz(-1.9705551) q[2];
sx q[2];
rz(-2.2054889) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6650508) q[1];
sx q[1];
rz(-2.9908123) q[1];
sx q[1];
rz(1.3097179) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4753597) q[3];
sx q[3];
rz(-1.1403475) q[3];
sx q[3];
rz(-0.43678771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6123885) q[2];
sx q[2];
rz(-0.91823053) q[2];
sx q[2];
rz(-2.6017792) q[2];
rz(-2.9811033) q[3];
sx q[3];
rz(-1.5925946) q[3];
sx q[3];
rz(2.3314355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4983343) q[0];
sx q[0];
rz(-0.67688268) q[0];
sx q[0];
rz(0.13370378) q[0];
rz(1.6224104) q[1];
sx q[1];
rz(-1.2956053) q[1];
sx q[1];
rz(-2.349283) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80817079) q[0];
sx q[0];
rz(-1.9229888) q[0];
sx q[0];
rz(2.9461224) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4033264) q[2];
sx q[2];
rz(-2.9472137) q[2];
sx q[2];
rz(2.3718718) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.49323248) q[1];
sx q[1];
rz(-2.2244452) q[1];
sx q[1];
rz(2.3834989) q[1];
rz(-pi) q[2];
rz(-3.1411857) q[3];
sx q[3];
rz(-2.8414629) q[3];
sx q[3];
rz(-0.89891922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.85492674) q[2];
sx q[2];
rz(-2.4714405) q[2];
sx q[2];
rz(2.3112042) q[2];
rz(1.3487799) q[3];
sx q[3];
rz(-1.034779) q[3];
sx q[3];
rz(-0.59876284) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26375672) q[0];
sx q[0];
rz(-0.97665518) q[0];
sx q[0];
rz(2.7184955) q[0];
rz(-1.9844203) q[1];
sx q[1];
rz(-1.6559699) q[1];
sx q[1];
rz(2.5194936) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7579278) q[0];
sx q[0];
rz(-2.1970941) q[0];
sx q[0];
rz(-0.48547283) q[0];
x q[1];
rz(-2.9329002) q[2];
sx q[2];
rz(-1.9819753) q[2];
sx q[2];
rz(-1.5087939) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9378949) q[1];
sx q[1];
rz(-0.40784281) q[1];
sx q[1];
rz(1.6126584) q[1];
rz(-pi) q[2];
rz(2.29633) q[3];
sx q[3];
rz(-1.5912531) q[3];
sx q[3];
rz(1.4485698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.28079924) q[2];
sx q[2];
rz(-1.4192899) q[2];
sx q[2];
rz(2.5234176) q[2];
rz(-1.482796) q[3];
sx q[3];
rz(-1.7890472) q[3];
sx q[3];
rz(-1.6331204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9789155) q[0];
sx q[0];
rz(-1.3293043) q[0];
sx q[0];
rz(0.83258122) q[0];
rz(-2.816448) q[1];
sx q[1];
rz(-2.145576) q[1];
sx q[1];
rz(-3.0772193) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8125809) q[0];
sx q[0];
rz(-1.6906749) q[0];
sx q[0];
rz(-2.8146312) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60026786) q[2];
sx q[2];
rz(-1.1986599) q[2];
sx q[2];
rz(1.063571) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1015013) q[1];
sx q[1];
rz(-1.027193) q[1];
sx q[1];
rz(0.60740791) q[1];
x q[2];
rz(-2.1983002) q[3];
sx q[3];
rz(-1.8648913) q[3];
sx q[3];
rz(2.7315745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9084106) q[2];
sx q[2];
rz(-2.5711214) q[2];
sx q[2];
rz(1.849966) q[2];
rz(2.6455247) q[3];
sx q[3];
rz(-0.78947624) q[3];
sx q[3];
rz(-1.1424278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.9555776) q[0];
sx q[0];
rz(-2.8514974) q[0];
sx q[0];
rz(0.41906038) q[0];
rz(-0.31173197) q[1];
sx q[1];
rz(-1.5429976) q[1];
sx q[1];
rz(-2.9980803) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1868033) q[0];
sx q[0];
rz(-1.1214897) q[0];
sx q[0];
rz(-1.9687998) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3355876) q[2];
sx q[2];
rz(-1.303788) q[2];
sx q[2];
rz(2.5036204) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4048481) q[1];
sx q[1];
rz(-2.2962034) q[1];
sx q[1];
rz(1.1400229) q[1];
rz(-1.7147948) q[3];
sx q[3];
rz(-1.8282951) q[3];
sx q[3];
rz(0.3443623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.781337) q[2];
sx q[2];
rz(-0.90136734) q[2];
sx q[2];
rz(1.1330913) q[2];
rz(-1.1059149) q[3];
sx q[3];
rz(-1.4191041) q[3];
sx q[3];
rz(-0.47479409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.356242) q[0];
sx q[0];
rz(-0.79371912) q[0];
sx q[0];
rz(1.445795) q[0];
rz(-0.71714199) q[1];
sx q[1];
rz(-1.8627867) q[1];
sx q[1];
rz(0.76146567) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26869795) q[0];
sx q[0];
rz(-2.1958662) q[0];
sx q[0];
rz(-1.752466) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0324208) q[2];
sx q[2];
rz(-2.3512977) q[2];
sx q[2];
rz(0.99144713) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.58757979) q[1];
sx q[1];
rz(-1.5214458) q[1];
sx q[1];
rz(-2.1448936) q[1];
rz(2.981212) q[3];
sx q[3];
rz(-1.3850936) q[3];
sx q[3];
rz(1.4726313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7722499) q[2];
sx q[2];
rz(-2.6476634) q[2];
sx q[2];
rz(-0.93144766) q[2];
rz(-2.5233968) q[3];
sx q[3];
rz(-1.5074916) q[3];
sx q[3];
rz(2.1759694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6922927) q[0];
sx q[0];
rz(-1.7789142) q[0];
sx q[0];
rz(-1.8071254) q[0];
rz(-2.6284699) q[1];
sx q[1];
rz(-0.98315364) q[1];
sx q[1];
rz(-1.2293053) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2199101) q[0];
sx q[0];
rz(-1.9905042) q[0];
sx q[0];
rz(2.3786663) q[0];
x q[1];
rz(-1.9077577) q[2];
sx q[2];
rz(-0.2744199) q[2];
sx q[2];
rz(-1.929259) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.25101623) q[1];
sx q[1];
rz(-0.15726798) q[1];
sx q[1];
rz(0.1319488) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5531814) q[3];
sx q[3];
rz(-2.2329139) q[3];
sx q[3];
rz(-2.3603242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8339771) q[2];
sx q[2];
rz(-0.5270842) q[2];
sx q[2];
rz(-2.9311467) q[2];
rz(3.0757507) q[3];
sx q[3];
rz(-2.6688771) q[3];
sx q[3];
rz(2.3426447) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54302067) q[0];
sx q[0];
rz(-0.6544756) q[0];
sx q[0];
rz(1.8684813) q[0];
rz(1.2741362) q[1];
sx q[1];
rz(-2.4576063) q[1];
sx q[1];
rz(-0.88409105) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9724204) q[0];
sx q[0];
rz(-2.5247249) q[0];
sx q[0];
rz(-1.5712117) q[0];
rz(-2.7124518) q[2];
sx q[2];
rz(-0.56641266) q[2];
sx q[2];
rz(2.8265116) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6581304) q[1];
sx q[1];
rz(-0.50848714) q[1];
sx q[1];
rz(2.1668651) q[1];
rz(-pi) q[2];
rz(-1.3119013) q[3];
sx q[3];
rz(-1.2719385) q[3];
sx q[3];
rz(-1.3204624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1511496) q[2];
sx q[2];
rz(-2.1572025) q[2];
sx q[2];
rz(-1.7669862) q[2];
rz(2.9511342) q[3];
sx q[3];
rz(-0.74646598) q[3];
sx q[3];
rz(-1.1838574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16354887) q[0];
sx q[0];
rz(-1.4530285) q[0];
sx q[0];
rz(1.2646041) q[0];
rz(-3.0679852) q[1];
sx q[1];
rz(-2.059748) q[1];
sx q[1];
rz(0.61990613) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093599565) q[0];
sx q[0];
rz(-1.7353829) q[0];
sx q[0];
rz(-2.3175146) q[0];
x q[1];
rz(-0.67155738) q[2];
sx q[2];
rz(-0.72973704) q[2];
sx q[2];
rz(-0.31873733) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.82295361) q[1];
sx q[1];
rz(-1.770449) q[1];
sx q[1];
rz(0.21491383) q[1];
x q[2];
rz(-1.8100753) q[3];
sx q[3];
rz(-1.115723) q[3];
sx q[3];
rz(-0.51067715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8668883) q[2];
sx q[2];
rz(-2.4394749) q[2];
sx q[2];
rz(0.14599027) q[2];
rz(0.22249666) q[3];
sx q[3];
rz(-2.9166418) q[3];
sx q[3];
rz(2.4116624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20919007) q[0];
sx q[0];
rz(-1.1840191) q[0];
sx q[0];
rz(-2.9970072) q[0];
rz(1.2296386) q[1];
sx q[1];
rz(-1.790779) q[1];
sx q[1];
rz(1.889224) q[1];
rz(1.7894921) q[2];
sx q[2];
rz(-0.98697294) q[2];
sx q[2];
rz(1.0728991) q[2];
rz(0.03224198) q[3];
sx q[3];
rz(-1.7809636) q[3];
sx q[3];
rz(-1.8973593) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
