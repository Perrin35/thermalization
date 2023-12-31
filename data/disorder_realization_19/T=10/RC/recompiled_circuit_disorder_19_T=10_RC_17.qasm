OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0269545) q[0];
sx q[0];
rz(-1.6898328) q[0];
sx q[0];
rz(0.64557689) q[0];
rz(-2.7627856) q[1];
sx q[1];
rz(-1.768755) q[1];
sx q[1];
rz(-1.6436613) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0092702) q[0];
sx q[0];
rz(-1.9267123) q[0];
sx q[0];
rz(0.15412553) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6662007) q[2];
sx q[2];
rz(-2.5089052) q[2];
sx q[2];
rz(-0.22104095) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9420535) q[1];
sx q[1];
rz(-0.76738165) q[1];
sx q[1];
rz(-2.8295423) q[1];
rz(-pi) q[2];
rz(1.9707768) q[3];
sx q[3];
rz(-2.1809289) q[3];
sx q[3];
rz(2.4657472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2215185) q[2];
sx q[2];
rz(-1.6211082) q[2];
sx q[2];
rz(-0.34040889) q[2];
rz(-0.83299625) q[3];
sx q[3];
rz(-0.70327988) q[3];
sx q[3];
rz(1.7849281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67265636) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(-0.54498589) q[0];
rz(-2.2333721) q[1];
sx q[1];
rz(-0.32749367) q[1];
sx q[1];
rz(-0.82495904) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23404113) q[0];
sx q[0];
rz(-1.3836765) q[0];
sx q[0];
rz(1.9812816) q[0];
rz(1.5812133) q[2];
sx q[2];
rz(-1.1679107) q[2];
sx q[2];
rz(-1.252623) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7488885) q[1];
sx q[1];
rz(-2.5671133) q[1];
sx q[1];
rz(-2.9267465) q[1];
rz(-2.1373848) q[3];
sx q[3];
rz(-1.7363747) q[3];
sx q[3];
rz(2.2752938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.58723441) q[2];
sx q[2];
rz(-2.6186269) q[2];
sx q[2];
rz(-2.9651802) q[2];
rz(-0.13088626) q[3];
sx q[3];
rz(-1.1754879) q[3];
sx q[3];
rz(3.1089354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.762887) q[0];
sx q[0];
rz(-0.58350199) q[0];
sx q[0];
rz(2.6718455) q[0];
rz(-1.6167971) q[1];
sx q[1];
rz(-0.47743118) q[1];
sx q[1];
rz(3.1030531) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3241203) q[0];
sx q[0];
rz(-1.5714374) q[0];
sx q[0];
rz(-1.5630432) q[0];
rz(-pi) q[1];
rz(2.5589587) q[2];
sx q[2];
rz(-2.8672672) q[2];
sx q[2];
rz(2.3290616) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1145648) q[1];
sx q[1];
rz(-1.6323043) q[1];
sx q[1];
rz(-2.0608749) q[1];
x q[2];
rz(2.3452957) q[3];
sx q[3];
rz(-2.3781804) q[3];
sx q[3];
rz(2.250716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5588351) q[2];
sx q[2];
rz(-1.0895412) q[2];
sx q[2];
rz(-2.2290686) q[2];
rz(-1.3085261) q[3];
sx q[3];
rz(-2.1378744) q[3];
sx q[3];
rz(-1.4413888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38347605) q[0];
sx q[0];
rz(-1.8261199) q[0];
sx q[0];
rz(0.34969774) q[0];
rz(1.2448467) q[1];
sx q[1];
rz(-0.30361509) q[1];
sx q[1];
rz(2.9464338) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53699025) q[0];
sx q[0];
rz(-1.450481) q[0];
sx q[0];
rz(-3.0483732) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85907016) q[2];
sx q[2];
rz(-0.67407437) q[2];
sx q[2];
rz(-0.40875834) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2745167) q[1];
sx q[1];
rz(-1.3356326) q[1];
sx q[1];
rz(-1.1067252) q[1];
x q[2];
rz(1.5870729) q[3];
sx q[3];
rz(-0.93173325) q[3];
sx q[3];
rz(1.0234327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23094709) q[2];
sx q[2];
rz(-0.4898943) q[2];
sx q[2];
rz(-1.267743) q[2];
rz(1.0686482) q[3];
sx q[3];
rz(-1.0125151) q[3];
sx q[3];
rz(-2.3012565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(1.1068263) q[0];
sx q[0];
rz(-1.6893457) q[0];
sx q[0];
rz(-0.1396133) q[0];
rz(1.0768249) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(0.37809125) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.970801) q[0];
sx q[0];
rz(-1.0266725) q[0];
sx q[0];
rz(2.8118954) q[0];
rz(0.81850448) q[2];
sx q[2];
rz(-2.2193925) q[2];
sx q[2];
rz(1.4865781) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2851583) q[1];
sx q[1];
rz(-0.57796961) q[1];
sx q[1];
rz(0.85698378) q[1];
rz(-1.7686754) q[3];
sx q[3];
rz(-1.1117001) q[3];
sx q[3];
rz(2.0898233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.87171626) q[2];
sx q[2];
rz(-1.3241974) q[2];
sx q[2];
rz(2.2536229) q[2];
rz(-0.97638431) q[3];
sx q[3];
rz(-1.7247) q[3];
sx q[3];
rz(-1.005727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3865005) q[0];
sx q[0];
rz(-3.0639102) q[0];
sx q[0];
rz(0.37762541) q[0];
rz(-0.32304421) q[1];
sx q[1];
rz(-2.4781365) q[1];
sx q[1];
rz(0.91517085) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5650425) q[0];
sx q[0];
rz(-2.8083907) q[0];
sx q[0];
rz(-2.8898426) q[0];
rz(-pi) q[1];
rz(-2.1377863) q[2];
sx q[2];
rz(-1.1319455) q[2];
sx q[2];
rz(-1.3958508) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9050161) q[1];
sx q[1];
rz(-1.6394776) q[1];
sx q[1];
rz(-2.828039) q[1];
x q[2];
rz(1.4573426) q[3];
sx q[3];
rz(-0.87943422) q[3];
sx q[3];
rz(-1.7951579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6043828) q[2];
sx q[2];
rz(-0.16468026) q[2];
sx q[2];
rz(-0.79052314) q[2];
rz(-2.8516155) q[3];
sx q[3];
rz(-2.4089549) q[3];
sx q[3];
rz(-2.524232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73615605) q[0];
sx q[0];
rz(-2.0744531) q[0];
sx q[0];
rz(0.36703584) q[0];
rz(1.908318) q[1];
sx q[1];
rz(-1.1739302) q[1];
sx q[1];
rz(-1.4253915) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4331626) q[0];
sx q[0];
rz(-1.3843952) q[0];
sx q[0];
rz(0.44048803) q[0];
rz(-pi) q[1];
rz(-3.0076214) q[2];
sx q[2];
rz(-2.6920762) q[2];
sx q[2];
rz(-1.9884895) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.094147625) q[1];
sx q[1];
rz(-0.74456753) q[1];
sx q[1];
rz(-2.5658957) q[1];
rz(-1.5951049) q[3];
sx q[3];
rz(-2.6377502) q[3];
sx q[3];
rz(-0.77038308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9795064) q[2];
sx q[2];
rz(-1.4922214) q[2];
sx q[2];
rz(-1.9308176) q[2];
rz(-3.1397505) q[3];
sx q[3];
rz(-0.76549923) q[3];
sx q[3];
rz(-2.4842998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2974671) q[0];
sx q[0];
rz(-2.5572889) q[0];
sx q[0];
rz(-3.0650744) q[0];
rz(2.466295) q[1];
sx q[1];
rz(-0.29390556) q[1];
sx q[1];
rz(0.75235596) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6589688) q[0];
sx q[0];
rz(-2.0730744) q[0];
sx q[0];
rz(1.9320022) q[0];
x q[1];
rz(0.92213995) q[2];
sx q[2];
rz(-1.4553918) q[2];
sx q[2];
rz(0.61308544) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.89531089) q[1];
sx q[1];
rz(-1.6591676) q[1];
sx q[1];
rz(1.6842711) q[1];
rz(-pi) q[2];
rz(-1.7767056) q[3];
sx q[3];
rz(-0.99109736) q[3];
sx q[3];
rz(1.0866144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.9383119) q[2];
sx q[2];
rz(-1.1802155) q[2];
sx q[2];
rz(-3.1414462) q[2];
rz(-1.1095307) q[3];
sx q[3];
rz(-0.90679449) q[3];
sx q[3];
rz(-2.8439567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9922441) q[0];
sx q[0];
rz(-0.23877564) q[0];
sx q[0];
rz(-1.0060271) q[0];
rz(-2.8736615) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(1.8267652) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6515503) q[0];
sx q[0];
rz(-2.393912) q[0];
sx q[0];
rz(2.7689395) q[0];
x q[1];
rz(2.8034231) q[2];
sx q[2];
rz(-0.38376946) q[2];
sx q[2];
rz(-1.9308117) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.73003522) q[1];
sx q[1];
rz(-0.81333232) q[1];
sx q[1];
rz(2.1557393) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3141765) q[3];
sx q[3];
rz(-1.5653059) q[3];
sx q[3];
rz(2.1237802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.43977794) q[2];
sx q[2];
rz(-0.54636991) q[2];
sx q[2];
rz(-0.24547274) q[2];
rz(-0.43073511) q[3];
sx q[3];
rz(-2.0499178) q[3];
sx q[3];
rz(0.47732863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4196639) q[0];
sx q[0];
rz(-0.99246445) q[0];
sx q[0];
rz(2.8549109) q[0];
rz(-0.87896705) q[1];
sx q[1];
rz(-1.6393839) q[1];
sx q[1];
rz(-0.39189664) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6673198) q[0];
sx q[0];
rz(-0.55756888) q[0];
sx q[0];
rz(-1.5865109) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91498615) q[2];
sx q[2];
rz(-1.6784385) q[2];
sx q[2];
rz(-2.8142625) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3956986) q[1];
sx q[1];
rz(-2.3682526) q[1];
sx q[1];
rz(1.8383154) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55962555) q[3];
sx q[3];
rz(-1.2979753) q[3];
sx q[3];
rz(0.19459693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.60951704) q[2];
sx q[2];
rz(-0.993002) q[2];
sx q[2];
rz(-1.9840476) q[2];
rz(-0.55084294) q[3];
sx q[3];
rz(-1.5778056) q[3];
sx q[3];
rz(1.1782066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75523238) q[0];
sx q[0];
rz(-1.7978783) q[0];
sx q[0];
rz(-1.88301) q[0];
rz(-1.8023087) q[1];
sx q[1];
rz(-2.5279999) q[1];
sx q[1];
rz(-2.7816714) q[1];
rz(2.9677283) q[2];
sx q[2];
rz(-2.3542913) q[2];
sx q[2];
rz(-2.3402294) q[2];
rz(-2.7950381) q[3];
sx q[3];
rz(-1.0481491) q[3];
sx q[3];
rz(-0.92878503) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
