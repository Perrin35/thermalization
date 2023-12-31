OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9392202) q[0];
sx q[0];
rz(-0.4063172) q[0];
sx q[0];
rz(-2.321474) q[0];
rz(2.7804873) q[1];
sx q[1];
rz(-0.63280025) q[1];
sx q[1];
rz(-0.83067218) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2031189) q[0];
sx q[0];
rz(-1.1955368) q[0];
sx q[0];
rz(-1.24036) q[0];
x q[1];
rz(2.6248706) q[2];
sx q[2];
rz(-0.80846723) q[2];
sx q[2];
rz(2.7441623) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9484529) q[1];
sx q[1];
rz(-1.2558736) q[1];
sx q[1];
rz(-0.83955168) q[1];
rz(-pi) q[2];
rz(2.6942263) q[3];
sx q[3];
rz(-1.9379741) q[3];
sx q[3];
rz(2.6255053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7044907) q[2];
sx q[2];
rz(-0.4318684) q[2];
sx q[2];
rz(3.0214156) q[2];
rz(-1.9834571) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(-2.2226298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0607818) q[0];
sx q[0];
rz(-1.7827543) q[0];
sx q[0];
rz(0.91180116) q[0];
rz(-0.78951019) q[1];
sx q[1];
rz(-0.98840886) q[1];
sx q[1];
rz(-2.8149014) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17900285) q[0];
sx q[0];
rz(-1.0674745) q[0];
sx q[0];
rz(-0.83053629) q[0];
x q[1];
rz(-0.15906449) q[2];
sx q[2];
rz(-2.0098364) q[2];
sx q[2];
rz(-2.1829) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.28194004) q[1];
sx q[1];
rz(-1.9743866) q[1];
sx q[1];
rz(2.2380026) q[1];
rz(3.0558415) q[3];
sx q[3];
rz(-0.88480703) q[3];
sx q[3];
rz(-0.12967295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.5102753) q[2];
sx q[2];
rz(-2.3687506) q[2];
sx q[2];
rz(1.2871683) q[2];
rz(0.10989799) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(1.7597642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.592955) q[0];
sx q[0];
rz(-0.73919636) q[0];
sx q[0];
rz(-0.32989311) q[0];
rz(-0.27711162) q[1];
sx q[1];
rz(-1.3169293) q[1];
sx q[1];
rz(2.0842016) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6431417) q[0];
sx q[0];
rz(-2.8006449) q[0];
sx q[0];
rz(1.5091512) q[0];
rz(-pi) q[1];
rz(-1.7705275) q[2];
sx q[2];
rz(-1.235851) q[2];
sx q[2];
rz(-1.0945601) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.33144618) q[1];
sx q[1];
rz(-2.3815037) q[1];
sx q[1];
rz(2.0600832) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3789165) q[3];
sx q[3];
rz(-1.9865611) q[3];
sx q[3];
rz(-0.92431812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.7827591) q[2];
sx q[2];
rz(-2.0262572) q[2];
sx q[2];
rz(-1.3624297) q[2];
rz(-0.6247012) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(-2.3220298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3574922) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(2.6065361) q[0];
rz(-2.0013981) q[1];
sx q[1];
rz(-1.813872) q[1];
sx q[1];
rz(-0.16539703) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9388158) q[0];
sx q[0];
rz(-1.446615) q[0];
sx q[0];
rz(-1.3846272) q[0];
rz(-pi) q[1];
rz(2.9781614) q[2];
sx q[2];
rz(-1.7497517) q[2];
sx q[2];
rz(2.2043259) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.71388984) q[1];
sx q[1];
rz(-1.0298924) q[1];
sx q[1];
rz(-1.9546024) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9966647) q[3];
sx q[3];
rz(-2.5440359) q[3];
sx q[3];
rz(-1.715341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7455204) q[2];
sx q[2];
rz(-1.3362276) q[2];
sx q[2];
rz(-0.20544927) q[2];
rz(1.127634) q[3];
sx q[3];
rz(-1.1536359) q[3];
sx q[3];
rz(1.8849461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(0.66185343) q[0];
sx q[0];
rz(-2.2275708) q[0];
sx q[0];
rz(-0.35476312) q[0];
rz(1.154249) q[1];
sx q[1];
rz(-2.216979) q[1];
sx q[1];
rz(-0.23194557) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97868279) q[0];
sx q[0];
rz(-1.9918348) q[0];
sx q[0];
rz(2.1091503) q[0];
rz(-pi) q[1];
rz(2.4460692) q[2];
sx q[2];
rz(-1.4624274) q[2];
sx q[2];
rz(2.7930789) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.31084727) q[1];
sx q[1];
rz(-1.7880882) q[1];
sx q[1];
rz(-0.22549916) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7730745) q[3];
sx q[3];
rz(-1.5526062) q[3];
sx q[3];
rz(0.25572488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0014687) q[2];
sx q[2];
rz(-1.3342369) q[2];
sx q[2];
rz(-1.2403963) q[2];
rz(0.59605789) q[3];
sx q[3];
rz(-1.8362935) q[3];
sx q[3];
rz(-1.8381455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0971138) q[0];
sx q[0];
rz(-0.070274027) q[0];
sx q[0];
rz(2.9388359) q[0];
rz(-0.98908201) q[1];
sx q[1];
rz(-1.443807) q[1];
sx q[1];
rz(0.99745497) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049012262) q[0];
sx q[0];
rz(-2.2708587) q[0];
sx q[0];
rz(2.9655365) q[0];
rz(2.7315797) q[2];
sx q[2];
rz(-2.3668681) q[2];
sx q[2];
rz(-2.0395525) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3012078) q[1];
sx q[1];
rz(-1.9921229) q[1];
sx q[1];
rz(-2.8788484) q[1];
rz(1.0783844) q[3];
sx q[3];
rz(-0.78740722) q[3];
sx q[3];
rz(-1.2315962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.22770195) q[2];
sx q[2];
rz(-1.9753549) q[2];
sx q[2];
rz(0.4513936) q[2];
rz(0.40870062) q[3];
sx q[3];
rz(-1.5390076) q[3];
sx q[3];
rz(-1.9394978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.8354427) q[0];
sx q[0];
rz(-2.7375484) q[0];
sx q[0];
rz(2.5174482) q[0];
rz(1.5165326) q[1];
sx q[1];
rz(-2.8846278) q[1];
sx q[1];
rz(0.61378941) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3261624) q[0];
sx q[0];
rz(-2.3809732) q[0];
sx q[0];
rz(-1.2994231) q[0];
x q[1];
rz(-2.9250547) q[2];
sx q[2];
rz(-2.8515186) q[2];
sx q[2];
rz(1.1494344) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.504096) q[1];
sx q[1];
rz(-1.4667965) q[1];
sx q[1];
rz(-1.6297479) q[1];
rz(-pi) q[2];
rz(1.8011439) q[3];
sx q[3];
rz(-0.47603546) q[3];
sx q[3];
rz(-1.3766833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7081786) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(1.1748574) q[2];
rz(-0.60837778) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(-1.7181989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90010086) q[0];
sx q[0];
rz(-1.8429723) q[0];
sx q[0];
rz(-0.4883782) q[0];
rz(1.6237367) q[1];
sx q[1];
rz(-1.3987712) q[1];
sx q[1];
rz(-0.98446313) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7229066) q[0];
sx q[0];
rz(-1.5333999) q[0];
sx q[0];
rz(-2.5953369) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1659053) q[2];
sx q[2];
rz(-2.3003909) q[2];
sx q[2];
rz(-1.0798432) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.24482803) q[1];
sx q[1];
rz(-2.7762189) q[1];
sx q[1];
rz(3.1246964) q[1];
x q[2];
rz(1.2582614) q[3];
sx q[3];
rz(-1.1141889) q[3];
sx q[3];
rz(1.6707591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.70665923) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(2.6064176) q[2];
rz(-2.0914071) q[3];
sx q[3];
rz(-1.340056) q[3];
sx q[3];
rz(-1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.500279) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(0.6859268) q[0];
rz(2.7507239) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(2.2156782) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5110916) q[0];
sx q[0];
rz(-1.7112268) q[0];
sx q[0];
rz(-1.3569843) q[0];
x q[1];
rz(-2.2885867) q[2];
sx q[2];
rz(-2.2249613) q[2];
sx q[2];
rz(0.6086364) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1123062) q[1];
sx q[1];
rz(-2.6767113) q[1];
sx q[1];
rz(-1.9025365) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20415281) q[3];
sx q[3];
rz(-1.6870058) q[3];
sx q[3];
rz(1.0128563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.19568504) q[2];
sx q[2];
rz(-0.47161272) q[2];
sx q[2];
rz(0.53608981) q[2];
rz(2.7219971) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(1.4887811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.0855899) q[0];
sx q[0];
rz(-0.35878006) q[0];
sx q[0];
rz(-2.7465903) q[0];
rz(-1.5123873) q[1];
sx q[1];
rz(-1.2724266) q[1];
sx q[1];
rz(-1.013247) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2241867) q[0];
sx q[0];
rz(-1.7099483) q[0];
sx q[0];
rz(1.805472) q[0];
x q[1];
rz(-0.11864885) q[2];
sx q[2];
rz(-0.35602202) q[2];
sx q[2];
rz(-0.3686541) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3328247) q[1];
sx q[1];
rz(-1.1507532) q[1];
sx q[1];
rz(-1.4648449) q[1];
rz(0.89684422) q[3];
sx q[3];
rz(-0.92771155) q[3];
sx q[3];
rz(0.53900063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44832486) q[2];
sx q[2];
rz(-1.6222745) q[2];
sx q[2];
rz(0.75941336) q[2];
rz(1.7761207) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(2.571648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41291819) q[0];
sx q[0];
rz(-1.2644132) q[0];
sx q[0];
rz(1.2046474) q[0];
rz(0.57327523) q[1];
sx q[1];
rz(-1.9075867) q[1];
sx q[1];
rz(1.8214068) q[1];
rz(0.076889597) q[2];
sx q[2];
rz(-2.0145609) q[2];
sx q[2];
rz(-1.2876074) q[2];
rz(-3.0782386) q[3];
sx q[3];
rz(-0.8739211) q[3];
sx q[3];
rz(-2.7779761) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
