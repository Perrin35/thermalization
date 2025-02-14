OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4545257) q[0];
sx q[0];
rz(5.0205686) q[0];
sx q[0];
rz(10.159259) q[0];
rz(1.6425411) q[1];
sx q[1];
rz(1.8416815) q[1];
sx q[1];
rz(11.587974) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9578732) q[0];
sx q[0];
rz(-1.5786202) q[0];
sx q[0];
rz(1.4343626) q[0];
x q[1];
rz(2.4327203) q[2];
sx q[2];
rz(-1.9908496) q[2];
sx q[2];
rz(1.5576897) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1803127) q[1];
sx q[1];
rz(-0.35113564) q[1];
sx q[1];
rz(-1.9957955) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4250578) q[3];
sx q[3];
rz(-0.50299931) q[3];
sx q[3];
rz(2.3921086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3413099) q[2];
sx q[2];
rz(-1.492123) q[2];
sx q[2];
rz(2.4260803) q[2];
rz(1.7320775) q[3];
sx q[3];
rz(-2.2655497) q[3];
sx q[3];
rz(1.5040262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8228804) q[0];
sx q[0];
rz(-0.57682288) q[0];
sx q[0];
rz(-3.1067644) q[0];
rz(2.4233129) q[1];
sx q[1];
rz(-1.4052582) q[1];
sx q[1];
rz(-1.1911596) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9485735) q[0];
sx q[0];
rz(-2.5002242) q[0];
sx q[0];
rz(1.8383584) q[0];
x q[1];
rz(-2.9600675) q[2];
sx q[2];
rz(-0.6904389) q[2];
sx q[2];
rz(-1.7275563) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.15219608) q[1];
sx q[1];
rz(-2.4862444) q[1];
sx q[1];
rz(2.6776777) q[1];
rz(-pi) q[2];
rz(0.19050772) q[3];
sx q[3];
rz(-1.1393875) q[3];
sx q[3];
rz(1.0781168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8027773) q[2];
sx q[2];
rz(-1.2884527) q[2];
sx q[2];
rz(3.0954933) q[2];
rz(2.3378546) q[3];
sx q[3];
rz(-1.9019144) q[3];
sx q[3];
rz(-2.5900335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4285202) q[0];
sx q[0];
rz(-1.5621194) q[0];
sx q[0];
rz(-0.61306104) q[0];
rz(-2.8763981) q[1];
sx q[1];
rz(-1.801633) q[1];
sx q[1];
rz(-0.048096098) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0408024) q[0];
sx q[0];
rz(-1.5682966) q[0];
sx q[0];
rz(2.0526171) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5597495) q[2];
sx q[2];
rz(-1.1238691) q[2];
sx q[2];
rz(2.4012411) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.887943) q[1];
sx q[1];
rz(-3.0620975) q[1];
sx q[1];
rz(-2.2924533) q[1];
x q[2];
rz(0.35086029) q[3];
sx q[3];
rz(-1.8985629) q[3];
sx q[3];
rz(1.3223437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4464104) q[2];
sx q[2];
rz(-1.1626652) q[2];
sx q[2];
rz(0.41333684) q[2];
rz(0.6684331) q[3];
sx q[3];
rz(-1.3276525) q[3];
sx q[3];
rz(1.7641164) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17871019) q[0];
sx q[0];
rz(-2.9229735) q[0];
sx q[0];
rz(0.14337732) q[0];
rz(-0.42477056) q[1];
sx q[1];
rz(-1.2062585) q[1];
sx q[1];
rz(-2.7075148) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.804337) q[0];
sx q[0];
rz(-1.4775663) q[0];
sx q[0];
rz(-3.1033959) q[0];
x q[1];
rz(-1.6878032) q[2];
sx q[2];
rz(-2.2196349) q[2];
sx q[2];
rz(-1.1800571) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0632759) q[1];
sx q[1];
rz(-0.5393636) q[1];
sx q[1];
rz(1.2961948) q[1];
rz(-pi) q[2];
rz(-0.8114154) q[3];
sx q[3];
rz(-1.7219117) q[3];
sx q[3];
rz(-0.18157696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.268959) q[2];
sx q[2];
rz(-0.7577529) q[2];
sx q[2];
rz(-0.41716519) q[2];
rz(-0.68552351) q[3];
sx q[3];
rz(-1.439582) q[3];
sx q[3];
rz(3.0912257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.3454469) q[0];
sx q[0];
rz(-0.31679994) q[0];
sx q[0];
rz(-1.5078804) q[0];
rz(-0.66868526) q[1];
sx q[1];
rz(-1.94328) q[1];
sx q[1];
rz(-2.0100458) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.749866) q[0];
sx q[0];
rz(-1.5380812) q[0];
sx q[0];
rz(-0.97264887) q[0];
rz(-pi) q[1];
rz(2.0756007) q[2];
sx q[2];
rz(-2.5194296) q[2];
sx q[2];
rz(1.0823859) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.9100489) q[1];
sx q[1];
rz(-2.0374415) q[1];
sx q[1];
rz(-1.779516) q[1];
x q[2];
rz(-0.82238166) q[3];
sx q[3];
rz(-0.17381343) q[3];
sx q[3];
rz(1.2863152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9097627) q[2];
sx q[2];
rz(-2.6015687) q[2];
sx q[2];
rz(-0.46169272) q[2];
rz(2.8716904) q[3];
sx q[3];
rz(-1.2609127) q[3];
sx q[3];
rz(0.65129483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57399026) q[0];
sx q[0];
rz(-0.42581588) q[0];
sx q[0];
rz(2.0800166) q[0];
rz(2.4827982) q[1];
sx q[1];
rz(-1.7196722) q[1];
sx q[1];
rz(0.40491358) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8773738) q[0];
sx q[0];
rz(-1.4658017) q[0];
sx q[0];
rz(-0.70606249) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54203029) q[2];
sx q[2];
rz(-1.2093636) q[2];
sx q[2];
rz(-2.3724477) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8078708) q[1];
sx q[1];
rz(-2.0757339) q[1];
sx q[1];
rz(-2.9823279) q[1];
x q[2];
rz(-0.93308927) q[3];
sx q[3];
rz(-1.328506) q[3];
sx q[3];
rz(-1.0090431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0460661) q[2];
sx q[2];
rz(-1.0045071) q[2];
sx q[2];
rz(-1.3859762) q[2];
rz(-2.4447377) q[3];
sx q[3];
rz(-2.160852) q[3];
sx q[3];
rz(2.3314886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54414576) q[0];
sx q[0];
rz(-0.61328855) q[0];
sx q[0];
rz(2.5469653) q[0];
rz(-1.1320629) q[1];
sx q[1];
rz(-1.4040399) q[1];
sx q[1];
rz(-0.51116699) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0027255) q[0];
sx q[0];
rz(-0.99926939) q[0];
sx q[0];
rz(-1.5062499) q[0];
x q[1];
rz(-2.2726502) q[2];
sx q[2];
rz(-0.41837439) q[2];
sx q[2];
rz(1.6623896) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8193002) q[1];
sx q[1];
rz(-0.27931133) q[1];
sx q[1];
rz(-1.929024) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53494142) q[3];
sx q[3];
rz(-2.3341093) q[3];
sx q[3];
rz(0.51715467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1749997) q[2];
sx q[2];
rz(-0.24093691) q[2];
sx q[2];
rz(-2.8540376) q[2];
rz(-1.1047085) q[3];
sx q[3];
rz(-1.46773) q[3];
sx q[3];
rz(-3.0159359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48968807) q[0];
sx q[0];
rz(-0.37864417) q[0];
sx q[0];
rz(2.4878159) q[0];
rz(-0.50103029) q[1];
sx q[1];
rz(-0.72622314) q[1];
sx q[1];
rz(2.3562866) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8619072) q[0];
sx q[0];
rz(-0.58478776) q[0];
sx q[0];
rz(-0.25081046) q[0];
x q[1];
rz(0.64898934) q[2];
sx q[2];
rz(-0.67399281) q[2];
sx q[2];
rz(1.6507932) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0846906) q[1];
sx q[1];
rz(-0.71699079) q[1];
sx q[1];
rz(0.3221953) q[1];
rz(-pi) q[2];
rz(-0.33881163) q[3];
sx q[3];
rz(-2.8991845) q[3];
sx q[3];
rz(0.95021866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.952897) q[2];
sx q[2];
rz(-1.6152103) q[2];
sx q[2];
rz(-0.3696028) q[2];
rz(-0.26436198) q[3];
sx q[3];
rz(-1.8601469) q[3];
sx q[3];
rz(3.1407691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0109176) q[0];
sx q[0];
rz(-0.70931804) q[0];
sx q[0];
rz(1.5189019) q[0];
rz(1.9021289) q[1];
sx q[1];
rz(-1.112097) q[1];
sx q[1];
rz(-2.1942203) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1110052) q[0];
sx q[0];
rz(-1.1717142) q[0];
sx q[0];
rz(-0.066094769) q[0];
x q[1];
rz(-2.9765509) q[2];
sx q[2];
rz(-2.7758935) q[2];
sx q[2];
rz(2.5132881) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4650146) q[1];
sx q[1];
rz(-1.7316707) q[1];
sx q[1];
rz(1.4032928) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2227433) q[3];
sx q[3];
rz(-0.17379856) q[3];
sx q[3];
rz(-0.17979515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6685247) q[2];
sx q[2];
rz(-0.99978414) q[2];
sx q[2];
rz(2.0295985) q[2];
rz(0.2019349) q[3];
sx q[3];
rz(-2.5578121) q[3];
sx q[3];
rz(-2.7659265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19875459) q[0];
sx q[0];
rz(-2.9534979) q[0];
sx q[0];
rz(-1.6756219) q[0];
rz(-1.6000043) q[1];
sx q[1];
rz(-1.8406638) q[1];
sx q[1];
rz(0.25513729) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45642325) q[0];
sx q[0];
rz(-3.0067634) q[0];
sx q[0];
rz(-1.441787) q[0];
rz(-pi) q[1];
rz(0.33619426) q[2];
sx q[2];
rz(-2.3393235) q[2];
sx q[2];
rz(0.65450571) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.49225351) q[1];
sx q[1];
rz(-1.6494177) q[1];
sx q[1];
rz(-1.2921974) q[1];
rz(-pi) q[2];
rz(-0.76948036) q[3];
sx q[3];
rz(-1.4325241) q[3];
sx q[3];
rz(1.6618123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.61350358) q[2];
sx q[2];
rz(-1.1897503) q[2];
sx q[2];
rz(2.7733754) q[2];
rz(-2.8585989) q[3];
sx q[3];
rz(-2.8734983) q[3];
sx q[3];
rz(-0.32576573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1267283) q[0];
sx q[0];
rz(-2.2751502) q[0];
sx q[0];
rz(-1.0979102) q[0];
rz(-0.49900838) q[1];
sx q[1];
rz(-1.9022763) q[1];
sx q[1];
rz(0.3872445) q[1];
rz(-1.6989742) q[2];
sx q[2];
rz(-1.2405996) q[2];
sx q[2];
rz(2.2157359) q[2];
rz(-1.3697237) q[3];
sx q[3];
rz(-1.5101505) q[3];
sx q[3];
rz(2.8240614) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
