OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20237246) q[0];
sx q[0];
rz(-2.7352754) q[0];
sx q[0];
rz(2.321474) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(0.63280025) q[1];
sx q[1];
rz(11.735698) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18618628) q[0];
sx q[0];
rz(-2.6468228) q[0];
sx q[0];
rz(-0.68899378) q[0];
rz(-0.73859282) q[2];
sx q[2];
rz(-1.2054218) q[2];
sx q[2];
rz(1.5473168) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.1931397) q[1];
sx q[1];
rz(-1.8857191) q[1];
sx q[1];
rz(-0.83955168) q[1];
x q[2];
rz(-1.9740231) q[3];
sx q[3];
rz(-1.1551757) q[3];
sx q[3];
rz(0.88413873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7044907) q[2];
sx q[2];
rz(-0.4318684) q[2];
sx q[2];
rz(-3.0214156) q[2];
rz(-1.9834571) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(0.91896287) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0607818) q[0];
sx q[0];
rz(-1.7827543) q[0];
sx q[0];
rz(2.2297915) q[0];
rz(-2.3520825) q[1];
sx q[1];
rz(-2.1531838) q[1];
sx q[1];
rz(-2.8149014) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90601901) q[0];
sx q[0];
rz(-0.86750194) q[0];
sx q[0];
rz(0.88615449) q[0];
x q[1];
rz(2.9825282) q[2];
sx q[2];
rz(-2.0098364) q[2];
sx q[2];
rz(-2.1829) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5527199) q[1];
sx q[1];
rz(-2.1761804) q[1];
sx q[1];
rz(2.6436716) q[1];
x q[2];
rz(0.88300206) q[3];
sx q[3];
rz(-1.5044754) q[3];
sx q[3];
rz(1.6460713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6313173) q[2];
sx q[2];
rz(-2.3687506) q[2];
sx q[2];
rz(-1.8544244) q[2];
rz(-0.10989799) q[3];
sx q[3];
rz(-1.7307614) q[3];
sx q[3];
rz(-1.7597642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(2.592955) q[0];
sx q[0];
rz(-0.73919636) q[0];
sx q[0];
rz(-2.8116995) q[0];
rz(0.27711162) q[1];
sx q[1];
rz(-1.8246633) q[1];
sx q[1];
rz(-1.057391) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6431417) q[0];
sx q[0];
rz(-0.34094778) q[0];
sx q[0];
rz(1.5091512) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3710652) q[2];
sx q[2];
rz(-1.9057416) q[2];
sx q[2];
rz(-2.0470326) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2705546) q[1];
sx q[1];
rz(-1.9005617) q[1];
sx q[1];
rz(2.2689181) q[1];
rz(-pi) q[2];
rz(1.7626761) q[3];
sx q[3];
rz(-1.9865611) q[3];
sx q[3];
rz(0.92431812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3588336) q[2];
sx q[2];
rz(-1.1153355) q[2];
sx q[2];
rz(1.7791629) q[2];
rz(-2.5168915) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(2.3220298) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7841004) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(0.5350565) q[0];
rz(-1.1401945) q[1];
sx q[1];
rz(-1.813872) q[1];
sx q[1];
rz(0.16539703) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9497313) q[0];
sx q[0];
rz(-2.9182069) q[0];
sx q[0];
rz(-2.1641157) q[0];
rz(-pi) q[1];
rz(2.303316) q[2];
sx q[2];
rz(-2.8998313) q[2];
sx q[2];
rz(-1.6844815) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.71388984) q[1];
sx q[1];
rz(-1.0298924) q[1];
sx q[1];
rz(-1.1869903) q[1];
x q[2];
rz(2.1256251) q[3];
sx q[3];
rz(-1.3362243) q[3];
sx q[3];
rz(-0.50333422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39607221) q[2];
sx q[2];
rz(-1.8053651) q[2];
sx q[2];
rz(0.20544927) q[2];
rz(1.127634) q[3];
sx q[3];
rz(-1.1536359) q[3];
sx q[3];
rz(-1.2566465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.4797392) q[0];
sx q[0];
rz(-0.91402188) q[0];
sx q[0];
rz(-2.7868295) q[0];
rz(1.154249) q[1];
sx q[1];
rz(-2.216979) q[1];
sx q[1];
rz(-0.23194557) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1922798) q[0];
sx q[0];
rz(-2.4711907) q[0];
sx q[0];
rz(0.85286661) q[0];
rz(-pi) q[1];
rz(-2.4460692) q[2];
sx q[2];
rz(-1.6791653) q[2];
sx q[2];
rz(-0.34851375) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9310589) q[1];
sx q[1];
rz(-1.3506883) q[1];
sx q[1];
rz(-1.7935497) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7730745) q[3];
sx q[3];
rz(-1.5889865) q[3];
sx q[3];
rz(2.8858678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.140124) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(-1.2403963) q[2];
rz(-0.59605789) q[3];
sx q[3];
rz(-1.3052992) q[3];
sx q[3];
rz(1.3034472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0444788) q[0];
sx q[0];
rz(-3.0713186) q[0];
sx q[0];
rz(2.9388359) q[0];
rz(2.1525106) q[1];
sx q[1];
rz(-1.443807) q[1];
sx q[1];
rz(-2.1441377) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7339242) q[0];
sx q[0];
rz(-1.4364388) q[0];
sx q[0];
rz(0.863048) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41001292) q[2];
sx q[2];
rz(-2.3668681) q[2];
sx q[2];
rz(-2.0395525) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.62085405) q[1];
sx q[1];
rz(-1.8100909) q[1];
sx q[1];
rz(-1.1362856) q[1];
x q[2];
rz(2.0632083) q[3];
sx q[3];
rz(-0.78740722) q[3];
sx q[3];
rz(1.2315962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.22770195) q[2];
sx q[2];
rz(-1.1662377) q[2];
sx q[2];
rz(2.690199) q[2];
rz(-0.40870062) q[3];
sx q[3];
rz(-1.6025851) q[3];
sx q[3];
rz(-1.9394978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8354427) q[0];
sx q[0];
rz(-0.4040443) q[0];
sx q[0];
rz(0.62414449) q[0];
rz(-1.5165326) q[1];
sx q[1];
rz(-2.8846278) q[1];
sx q[1];
rz(-0.61378941) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3261624) q[0];
sx q[0];
rz(-0.76061941) q[0];
sx q[0];
rz(-1.8421696) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8579312) q[2];
sx q[2];
rz(-1.5093056) q[2];
sx q[2];
rz(2.9279857) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.63749667) q[1];
sx q[1];
rz(-1.6747961) q[1];
sx q[1];
rz(-1.6297479) q[1];
rz(-3.0244163) q[3];
sx q[3];
rz(-1.1083318) q[3];
sx q[3];
rz(1.6346491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7081786) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(-1.1748574) q[2];
rz(0.60837778) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(-1.4233937) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90010086) q[0];
sx q[0];
rz(-1.8429723) q[0];
sx q[0];
rz(-0.4883782) q[0];
rz(-1.5178559) q[1];
sx q[1];
rz(-1.3987712) q[1];
sx q[1];
rz(2.1571295) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9280633) q[0];
sx q[0];
rz(-0.54740471) q[0];
sx q[0];
rz(0.071896032) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9756873) q[2];
sx q[2];
rz(-2.3003909) q[2];
sx q[2];
rz(2.0617495) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.24482803) q[1];
sx q[1];
rz(-0.36537376) q[1];
sx q[1];
rz(-3.1246964) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5823309) q[3];
sx q[3];
rz(-0.54702938) q[3];
sx q[3];
rz(0.83838851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4349334) q[2];
sx q[2];
rz(-1.6985396) q[2];
sx q[2];
rz(2.6064176) q[2];
rz(2.0914071) q[3];
sx q[3];
rz(-1.8015367) q[3];
sx q[3];
rz(-1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.500279) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(-2.4556659) q[0];
rz(0.39086875) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(-2.2156782) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5092963) q[0];
sx q[0];
rz(-0.25521454) q[0];
sx q[0];
rz(-0.98310982) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3472896) q[2];
sx q[2];
rz(-1.0216733) q[2];
sx q[2];
rz(1.6910451) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.33854252) q[1];
sx q[1];
rz(-2.0085137) q[1];
sx q[1];
rz(2.9796757) q[1];
x q[2];
rz(-2.6191606) q[3];
sx q[3];
rz(-2.9070832) q[3];
sx q[3];
rz(2.0731376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.19568504) q[2];
sx q[2];
rz(-0.47161272) q[2];
sx q[2];
rz(0.53608981) q[2];
rz(-0.4195956) q[3];
sx q[3];
rz(-2.5511238) q[3];
sx q[3];
rz(-1.6528116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0560028) q[0];
sx q[0];
rz(-2.7828126) q[0];
sx q[0];
rz(-0.39500239) q[0];
rz(1.6292054) q[1];
sx q[1];
rz(-1.869166) q[1];
sx q[1];
rz(1.013247) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9625044) q[0];
sx q[0];
rz(-0.27217406) q[0];
sx q[0];
rz(1.0286691) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6147862) q[2];
sx q[2];
rz(-1.2173875) q[2];
sx q[2];
rz(-0.24214889) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4229065) q[1];
sx q[1];
rz(-1.4740853) q[1];
sx q[1];
rz(0.42214091) q[1];
rz(0.89684422) q[3];
sx q[3];
rz(-2.2138811) q[3];
sx q[3];
rz(2.602592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6932678) q[2];
sx q[2];
rz(-1.5193181) q[2];
sx q[2];
rz(0.75941336) q[2];
rz(1.365472) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(0.56994462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7286745) q[0];
sx q[0];
rz(-1.8771794) q[0];
sx q[0];
rz(-1.9369453) q[0];
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
