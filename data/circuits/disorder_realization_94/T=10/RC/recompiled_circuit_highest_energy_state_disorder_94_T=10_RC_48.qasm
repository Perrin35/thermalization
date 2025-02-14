OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2430159) q[0];
sx q[0];
rz(-1.1216811) q[0];
sx q[0];
rz(0.46410528) q[0];
rz(2.3857181) q[1];
sx q[1];
rz(-0.90944374) q[1];
sx q[1];
rz(-1.3165201) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6073734) q[0];
sx q[0];
rz(-1.6053622) q[0];
sx q[0];
rz(0.065163356) q[0];
rz(1.6208956) q[2];
sx q[2];
rz(-1.230579) q[2];
sx q[2];
rz(0.17776793) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2006372) q[1];
sx q[1];
rz(-0.66087729) q[1];
sx q[1];
rz(2.6913068) q[1];
rz(-pi) q[2];
rz(0.25578292) q[3];
sx q[3];
rz(-1.2885416) q[3];
sx q[3];
rz(-1.5402067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.89062771) q[2];
sx q[2];
rz(-2.784412) q[2];
sx q[2];
rz(-0.29699057) q[2];
rz(2.4885528) q[3];
sx q[3];
rz(-1.5584471) q[3];
sx q[3];
rz(2.0853341) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62946573) q[0];
sx q[0];
rz(-2.1194206) q[0];
sx q[0];
rz(0.28489354) q[0];
rz(3.0902872) q[1];
sx q[1];
rz(-2.3785794) q[1];
sx q[1];
rz(2.5047393) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92572407) q[0];
sx q[0];
rz(-2.0179739) q[0];
sx q[0];
rz(2.607292) q[0];
rz(-pi) q[1];
rz(1.5696196) q[2];
sx q[2];
rz(-1.92244) q[2];
sx q[2];
rz(-0.52368173) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7015345) q[1];
sx q[1];
rz(-1.4949168) q[1];
sx q[1];
rz(-2.9666191) q[1];
x q[2];
rz(2.4356682) q[3];
sx q[3];
rz(-1.4691938) q[3];
sx q[3];
rz(1.3639579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9597943) q[2];
sx q[2];
rz(-0.85796285) q[2];
sx q[2];
rz(-2.174343) q[2];
rz(2.8773384) q[3];
sx q[3];
rz(-1.6382917) q[3];
sx q[3];
rz(-1.9056457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0018175) q[0];
sx q[0];
rz(-1.5622666) q[0];
sx q[0];
rz(-2.0215969) q[0];
rz(-2.6657875) q[1];
sx q[1];
rz(-2.0640524) q[1];
sx q[1];
rz(2.8376104) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33789269) q[0];
sx q[0];
rz(-0.56625444) q[0];
sx q[0];
rz(1.6223652) q[0];
x q[1];
rz(0.23786942) q[2];
sx q[2];
rz(-1.6697236) q[2];
sx q[2];
rz(1.7406657) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7870362) q[1];
sx q[1];
rz(-2.4156385) q[1];
sx q[1];
rz(-1.7257084) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.899053) q[3];
sx q[3];
rz(-0.69787301) q[3];
sx q[3];
rz(-2.1062811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1659282) q[2];
sx q[2];
rz(-2.1037481) q[2];
sx q[2];
rz(2.4456639) q[2];
rz(-1.1337229) q[3];
sx q[3];
rz(-1.1997831) q[3];
sx q[3];
rz(2.5886562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7021779) q[0];
sx q[0];
rz(-1.004383) q[0];
sx q[0];
rz(1.0473921) q[0];
rz(1.1737431) q[1];
sx q[1];
rz(-2.4672697) q[1];
sx q[1];
rz(-1.6204999) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3701909) q[0];
sx q[0];
rz(-2.6242497) q[0];
sx q[0];
rz(2.8964554) q[0];
x q[1];
rz(2.8212566) q[2];
sx q[2];
rz(-2.2514859) q[2];
sx q[2];
rz(2.5706511) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7340659) q[1];
sx q[1];
rz(-2.1868949) q[1];
sx q[1];
rz(-1.0259088) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32471809) q[3];
sx q[3];
rz(-1.8408662) q[3];
sx q[3];
rz(2.3809759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0461222) q[2];
sx q[2];
rz(-2.0746524) q[2];
sx q[2];
rz(-0.53786892) q[2];
rz(-1.6543903) q[3];
sx q[3];
rz(-2.7345149) q[3];
sx q[3];
rz(2.4552104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2122413) q[0];
sx q[0];
rz(-2.90726) q[0];
sx q[0];
rz(2.5382163) q[0];
rz(-2.3790908) q[1];
sx q[1];
rz(-1.1926032) q[1];
sx q[1];
rz(0.71294436) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.514587) q[0];
sx q[0];
rz(-2.7381185) q[0];
sx q[0];
rz(-1.7395942) q[0];
rz(-2.0335401) q[2];
sx q[2];
rz(-1.4456835) q[2];
sx q[2];
rz(2.3274328) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7984287) q[1];
sx q[1];
rz(-0.53939941) q[1];
sx q[1];
rz(-1.1718962) q[1];
x q[2];
rz(-0.74728031) q[3];
sx q[3];
rz(-0.76939121) q[3];
sx q[3];
rz(-1.4667433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.9731628) q[2];
sx q[2];
rz(-0.91732401) q[2];
sx q[2];
rz(-1.6707576) q[2];
rz(2.6245608) q[3];
sx q[3];
rz(-1.0443338) q[3];
sx q[3];
rz(1.8487336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9789199) q[0];
sx q[0];
rz(-0.51386583) q[0];
sx q[0];
rz(-0.55147076) q[0];
rz(3.1061213) q[1];
sx q[1];
rz(-2.008581) q[1];
sx q[1];
rz(-1.5303401) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0198675) q[0];
sx q[0];
rz(-2.1824153) q[0];
sx q[0];
rz(0.83329552) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6257203) q[2];
sx q[2];
rz(-2.3936205) q[2];
sx q[2];
rz(2.2750281) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4582461) q[1];
sx q[1];
rz(-2.2846232) q[1];
sx q[1];
rz(0.48269646) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2238268) q[3];
sx q[3];
rz(-1.5035331) q[3];
sx q[3];
rz(-1.3239087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0772721) q[2];
sx q[2];
rz(-2.6214226) q[2];
sx q[2];
rz(-1.8989829) q[2];
rz(2.6284435) q[3];
sx q[3];
rz(-0.41392252) q[3];
sx q[3];
rz(-1.7665524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4921017) q[0];
sx q[0];
rz(-2.212337) q[0];
sx q[0];
rz(-3.0249366) q[0];
rz(0.77313441) q[1];
sx q[1];
rz(-1.1583068) q[1];
sx q[1];
rz(-1.6366417) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7709383) q[0];
sx q[0];
rz(-2.2794834) q[0];
sx q[0];
rz(0.38972008) q[0];
rz(1.1820311) q[2];
sx q[2];
rz(-1.6259527) q[2];
sx q[2];
rz(-0.15832947) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4915858) q[1];
sx q[1];
rz(-1.9284231) q[1];
sx q[1];
rz(2.3050344) q[1];
x q[2];
rz(-0.36913334) q[3];
sx q[3];
rz(-0.68644409) q[3];
sx q[3];
rz(2.4001207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2152805) q[2];
sx q[2];
rz(-2.8947688) q[2];
sx q[2];
rz(-0.8872633) q[2];
rz(0.67982802) q[3];
sx q[3];
rz(-0.66893783) q[3];
sx q[3];
rz(-2.5086856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.585084) q[0];
sx q[0];
rz(-1.8483138) q[0];
sx q[0];
rz(0.65628091) q[0];
rz(2.9346924) q[1];
sx q[1];
rz(-1.921939) q[1];
sx q[1];
rz(-1.9237178) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042699634) q[0];
sx q[0];
rz(-1.4880565) q[0];
sx q[0];
rz(-1.6059884) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0711259) q[2];
sx q[2];
rz(-2.4917951) q[2];
sx q[2];
rz(-1.4184679) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7102573) q[1];
sx q[1];
rz(-1.6847622) q[1];
sx q[1];
rz(-0.15941391) q[1];
rz(-0.37533203) q[3];
sx q[3];
rz(-2.4960244) q[3];
sx q[3];
rz(-0.69824346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8021585) q[2];
sx q[2];
rz(-1.755244) q[2];
sx q[2];
rz(-2.459724) q[2];
rz(-1.1784461) q[3];
sx q[3];
rz(-2.384187) q[3];
sx q[3];
rz(1.9044378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7449529) q[0];
sx q[0];
rz(-1.6093901) q[0];
sx q[0];
rz(-2.9314281) q[0];
rz(-0.22178966) q[1];
sx q[1];
rz(-0.91786018) q[1];
sx q[1];
rz(-0.04235696) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6307459) q[0];
sx q[0];
rz(-1.5928942) q[0];
sx q[0];
rz(-2.9539897) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3826959) q[2];
sx q[2];
rz(-2.3504371) q[2];
sx q[2];
rz(3.0393785) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.90396229) q[1];
sx q[1];
rz(-2.1231066) q[1];
sx q[1];
rz(1.2321074) q[1];
rz(-pi) q[2];
rz(2.0324519) q[3];
sx q[3];
rz(-1.084514) q[3];
sx q[3];
rz(-1.2580296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.28045851) q[2];
sx q[2];
rz(-2.0529842) q[2];
sx q[2];
rz(-0.66656485) q[2];
rz(2.376453) q[3];
sx q[3];
rz(-2.9355526) q[3];
sx q[3];
rz(-1.7407181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6078981) q[0];
sx q[0];
rz(-2.7602637) q[0];
sx q[0];
rz(-0.68699849) q[0];
rz(-1.2835361) q[1];
sx q[1];
rz(-1.9891519) q[1];
sx q[1];
rz(-2.5680465) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4592996) q[0];
sx q[0];
rz(-2.1577765) q[0];
sx q[0];
rz(0.94265509) q[0];
rz(-2.0229265) q[2];
sx q[2];
rz(-1.7799449) q[2];
sx q[2];
rz(0.5224519) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.191205) q[1];
sx q[1];
rz(-2.5053686) q[1];
sx q[1];
rz(-2.1649182) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76400842) q[3];
sx q[3];
rz(-2.2413018) q[3];
sx q[3];
rz(1.7830199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0428697) q[2];
sx q[2];
rz(-1.1770153) q[2];
sx q[2];
rz(0.086611835) q[2];
rz(0.11836554) q[3];
sx q[3];
rz(-1.5425073) q[3];
sx q[3];
rz(-0.043244403) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31556986) q[0];
sx q[0];
rz(-2.12349) q[0];
sx q[0];
rz(0.77793599) q[0];
rz(-0.28620537) q[1];
sx q[1];
rz(-2.310391) q[1];
sx q[1];
rz(-1.8117767) q[1];
rz(2.8092842) q[2];
sx q[2];
rz(-1.7996017) q[2];
sx q[2];
rz(-1.9881291) q[2];
rz(1.8686915) q[3];
sx q[3];
rz(-1.1392698) q[3];
sx q[3];
rz(-1.4261693) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
