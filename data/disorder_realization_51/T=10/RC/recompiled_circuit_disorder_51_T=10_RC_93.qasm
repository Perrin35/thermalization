OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5470619) q[0];
sx q[0];
rz(4.2630258) q[0];
sx q[0];
rz(6.1017258) q[0];
rz(-1.0815066) q[1];
sx q[1];
rz(-2.4681611) q[1];
sx q[1];
rz(-2.0884617) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5232789) q[0];
sx q[0];
rz(-1.7821454) q[0];
sx q[0];
rz(2.0413415) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62280957) q[2];
sx q[2];
rz(-1.8138759) q[2];
sx q[2];
rz(-2.958975) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.92703544) q[1];
sx q[1];
rz(-1.7662449) q[1];
sx q[1];
rz(-0.77597015) q[1];
x q[2];
rz(2.0041204) q[3];
sx q[3];
rz(-2.7702799) q[3];
sx q[3];
rz(2.5759047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.16333214) q[2];
sx q[2];
rz(-2.0330727) q[2];
sx q[2];
rz(1.7738316) q[2];
rz(1.0129499) q[3];
sx q[3];
rz(-2.2949341) q[3];
sx q[3];
rz(-3.0701385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0435836) q[0];
sx q[0];
rz(-1.4494891) q[0];
sx q[0];
rz(-2.5464771) q[0];
rz(2.0960506) q[1];
sx q[1];
rz(-1.7273993) q[1];
sx q[1];
rz(-1.6275303) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1925416) q[0];
sx q[0];
rz(-2.4727614) q[0];
sx q[0];
rz(2.2632954) q[0];
rz(-pi) q[1];
rz(-0.83346955) q[2];
sx q[2];
rz(-1.5499299) q[2];
sx q[2];
rz(1.8881063) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0552057) q[1];
sx q[1];
rz(-1.7654164) q[1];
sx q[1];
rz(-0.17533949) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4684832) q[3];
sx q[3];
rz(-1.6204837) q[3];
sx q[3];
rz(-1.697584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.65511584) q[2];
sx q[2];
rz(-1.0412419) q[2];
sx q[2];
rz(-0.79616037) q[2];
rz(-2.1697309) q[3];
sx q[3];
rz(-0.704851) q[3];
sx q[3];
rz(0.17175737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3765091) q[0];
sx q[0];
rz(-2.3628545) q[0];
sx q[0];
rz(-3.0774975) q[0];
rz(-0.31072101) q[1];
sx q[1];
rz(-1.6711845) q[1];
sx q[1];
rz(1.4583189) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38760936) q[0];
sx q[0];
rz(-1.4235272) q[0];
sx q[0];
rz(-3.0940042) q[0];
rz(0.73297357) q[2];
sx q[2];
rz(-0.79967116) q[2];
sx q[2];
rz(-0.016499585) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.872936) q[1];
sx q[1];
rz(-1.5169414) q[1];
sx q[1];
rz(-0.40952803) q[1];
rz(2.8267616) q[3];
sx q[3];
rz(-1.587084) q[3];
sx q[3];
rz(-2.6967238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.98214275) q[2];
sx q[2];
rz(-0.84656707) q[2];
sx q[2];
rz(-0.67908755) q[2];
rz(-1.3726161) q[3];
sx q[3];
rz(-1.8434098) q[3];
sx q[3];
rz(0.94846559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8452334) q[0];
sx q[0];
rz(-0.56931749) q[0];
sx q[0];
rz(-2.0171719) q[0];
rz(1.2202948) q[1];
sx q[1];
rz(-1.1923469) q[1];
sx q[1];
rz(-1.6569998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35614466) q[0];
sx q[0];
rz(-1.6011642) q[0];
sx q[0];
rz(-1.6214451) q[0];
rz(2.0882656) q[2];
sx q[2];
rz(-1.7382858) q[2];
sx q[2];
rz(2.6094764) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2034519) q[1];
sx q[1];
rz(-2.1362937) q[1];
sx q[1];
rz(-2.6184665) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4184166) q[3];
sx q[3];
rz(-1.8225065) q[3];
sx q[3];
rz(2.6808006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1241887) q[2];
sx q[2];
rz(-1.2839395) q[2];
sx q[2];
rz(-1.0162639) q[2];
rz(1.7381564) q[3];
sx q[3];
rz(-1.6442464) q[3];
sx q[3];
rz(-2.0531634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6338585) q[0];
sx q[0];
rz(-1.2861179) q[0];
sx q[0];
rz(-2.741709) q[0];
rz(-1.1625066) q[1];
sx q[1];
rz(-1.8116654) q[1];
sx q[1];
rz(-0.17366017) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7090209) q[0];
sx q[0];
rz(-0.077843277) q[0];
sx q[0];
rz(-0.37065931) q[0];
rz(-pi) q[1];
rz(-1.7816254) q[2];
sx q[2];
rz(-1.1846917) q[2];
sx q[2];
rz(0.93036508) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9160737) q[1];
sx q[1];
rz(-1.5991296) q[1];
sx q[1];
rz(-1.4654935) q[1];
rz(-pi) q[2];
rz(1.7791622) q[3];
sx q[3];
rz(-1.8959798) q[3];
sx q[3];
rz(-1.9735826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6614723) q[2];
sx q[2];
rz(-1.9813333) q[2];
sx q[2];
rz(-2.373467) q[2];
rz(-0.85401946) q[3];
sx q[3];
rz(-1.7212399) q[3];
sx q[3];
rz(2.1693726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0356692) q[0];
sx q[0];
rz(-2.8978455) q[0];
sx q[0];
rz(1.4703898) q[0];
rz(-0.51180965) q[1];
sx q[1];
rz(-2.6302331) q[1];
sx q[1];
rz(1.9981729) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4758159) q[0];
sx q[0];
rz(-1.4205298) q[0];
sx q[0];
rz(-1.6554324) q[0];
rz(-2.5383699) q[2];
sx q[2];
rz(-1.8458741) q[2];
sx q[2];
rz(3.0657257) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9547564) q[1];
sx q[1];
rz(-1.3969159) q[1];
sx q[1];
rz(1.9435005) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8361736) q[3];
sx q[3];
rz(-1.050204) q[3];
sx q[3];
rz(-0.49314317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.64289552) q[2];
sx q[2];
rz(-0.84474793) q[2];
sx q[2];
rz(-1.6112304) q[2];
rz(1.6879843) q[3];
sx q[3];
rz(-2.0724824) q[3];
sx q[3];
rz(0.24916515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0221508) q[0];
sx q[0];
rz(-0.74437737) q[0];
sx q[0];
rz(1.4402333) q[0];
rz(0.72921905) q[1];
sx q[1];
rz(-1.9897285) q[1];
sx q[1];
rz(2.008332) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6183375) q[0];
sx q[0];
rz(-1.6132857) q[0];
sx q[0];
rz(0.32251127) q[0];
rz(-pi) q[1];
x q[1];
rz(0.015958162) q[2];
sx q[2];
rz(-2.4292813) q[2];
sx q[2];
rz(-2.6238837) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.28753528) q[1];
sx q[1];
rz(-2.4234728) q[1];
sx q[1];
rz(1.1487886) q[1];
x q[2];
rz(0.53448581) q[3];
sx q[3];
rz(-2.839698) q[3];
sx q[3];
rz(1.2045977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.40522727) q[2];
sx q[2];
rz(-1.1232802) q[2];
sx q[2];
rz(0.48842946) q[2];
rz(-1.3119665) q[3];
sx q[3];
rz(-3.0977111) q[3];
sx q[3];
rz(-0.64129889) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0768123) q[0];
sx q[0];
rz(-1.8490054) q[0];
sx q[0];
rz(3.0704165) q[0];
rz(3.1094303) q[1];
sx q[1];
rz(-1.8036489) q[1];
sx q[1];
rz(-1.9326928) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6541518) q[0];
sx q[0];
rz(-1.3432661) q[0];
sx q[0];
rz(2.3339416) q[0];
x q[1];
rz(-2.2279943) q[2];
sx q[2];
rz(-2.845394) q[2];
sx q[2];
rz(-1.2072472) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3370918) q[1];
sx q[1];
rz(-0.13131222) q[1];
sx q[1];
rz(-1.4485703) q[1];
rz(-2.6505359) q[3];
sx q[3];
rz(-0.68248442) q[3];
sx q[3];
rz(1.3852937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.20748392) q[2];
sx q[2];
rz(-2.948163) q[2];
sx q[2];
rz(-1.5709546) q[2];
rz(0.87336826) q[3];
sx q[3];
rz(-1.7375172) q[3];
sx q[3];
rz(2.7895555) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97380012) q[0];
sx q[0];
rz(-1.6163102) q[0];
sx q[0];
rz(-2.8299676) q[0];
rz(-0.82178003) q[1];
sx q[1];
rz(-2.5506134) q[1];
sx q[1];
rz(-1.5100381) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8845997) q[0];
sx q[0];
rz(-1.8879226) q[0];
sx q[0];
rz(1.7811799) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0636343) q[2];
sx q[2];
rz(-2.7535451) q[2];
sx q[2];
rz(-1.0840814) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1875302) q[1];
sx q[1];
rz(-0.64779753) q[1];
sx q[1];
rz(2.8026583) q[1];
rz(-0.23267965) q[3];
sx q[3];
rz(-1.3921229) q[3];
sx q[3];
rz(1.0089547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2733549) q[2];
sx q[2];
rz(-2.6028825) q[2];
sx q[2];
rz(-0.81595016) q[2];
rz(2.6319035) q[3];
sx q[3];
rz(-2.3563801) q[3];
sx q[3];
rz(-1.2379439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.8907392) q[0];
sx q[0];
rz(-1.4287404) q[0];
sx q[0];
rz(-2.9113286) q[0];
rz(-0.62581217) q[1];
sx q[1];
rz(-2.1964549) q[1];
sx q[1];
rz(-2.4831916) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.828016) q[0];
sx q[0];
rz(-0.99452924) q[0];
sx q[0];
rz(2.1616031) q[0];
rz(0.41215956) q[2];
sx q[2];
rz(-2.4240652) q[2];
sx q[2];
rz(1.621643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.60724466) q[1];
sx q[1];
rz(-1.4048647) q[1];
sx q[1];
rz(-2.0545309) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2213311) q[3];
sx q[3];
rz(-2.0798827) q[3];
sx q[3];
rz(-1.5290608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4663503) q[2];
sx q[2];
rz(-0.8224951) q[2];
sx q[2];
rz(-2.8038483) q[2];
rz(2.1181469) q[3];
sx q[3];
rz(-1.3600391) q[3];
sx q[3];
rz(-1.9256928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52453775) q[0];
sx q[0];
rz(-2.239997) q[0];
sx q[0];
rz(-0.064185113) q[0];
rz(1.2383923) q[1];
sx q[1];
rz(-1.7108142) q[1];
sx q[1];
rz(-1.6731813) q[1];
rz(2.9624883) q[2];
sx q[2];
rz(-2.1894107) q[2];
sx q[2];
rz(-1.4084963) q[2];
rz(3.0425439) q[3];
sx q[3];
rz(-2.5669813) q[3];
sx q[3];
rz(3.0909227) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];