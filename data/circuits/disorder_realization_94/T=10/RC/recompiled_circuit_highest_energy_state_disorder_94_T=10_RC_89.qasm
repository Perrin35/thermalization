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
rz(-2.6774874) q[0];
rz(2.3857181) q[1];
sx q[1];
rz(-0.90944374) q[1];
sx q[1];
rz(-1.3165201) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1027604) q[0];
sx q[0];
rz(-1.505672) q[0];
sx q[0];
rz(1.6054356) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0010536) q[2];
sx q[2];
rz(-0.34374434) q[2];
sx q[2];
rz(-0.02862169) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7500748) q[1];
sx q[1];
rz(-2.1563091) q[1];
sx q[1];
rz(1.2445009) q[1];
rz(-2.2881704) q[3];
sx q[3];
rz(-2.7629921) q[3];
sx q[3];
rz(-0.84747696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.89062771) q[2];
sx q[2];
rz(-0.35718063) q[2];
sx q[2];
rz(0.29699057) q[2];
rz(-2.4885528) q[3];
sx q[3];
rz(-1.5584471) q[3];
sx q[3];
rz(1.0562586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62946573) q[0];
sx q[0];
rz(-1.0221721) q[0];
sx q[0];
rz(2.8566991) q[0];
rz(3.0902872) q[1];
sx q[1];
rz(-2.3785794) q[1];
sx q[1];
rz(2.5047393) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2460265) q[0];
sx q[0];
rz(-1.0937397) q[0];
sx q[0];
rz(-1.0624079) q[0];
rz(0.35164386) q[2];
sx q[2];
rz(-1.571901) q[2];
sx q[2];
rz(1.0467093) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6057558) q[1];
sx q[1];
rz(-0.19056377) q[1];
sx q[1];
rz(0.41175731) q[1];
rz(-pi) q[2];
rz(-1.4376182) q[3];
sx q[3];
rz(-2.2723327) q[3];
sx q[3];
rz(-3.0210036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9597943) q[2];
sx q[2];
rz(-0.85796285) q[2];
sx q[2];
rz(0.96724969) q[2];
rz(0.26425427) q[3];
sx q[3];
rz(-1.6382917) q[3];
sx q[3];
rz(1.9056457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0018175) q[0];
sx q[0];
rz(-1.5793261) q[0];
sx q[0];
rz(-1.1199957) q[0];
rz(-2.6657875) q[1];
sx q[1];
rz(-1.0775403) q[1];
sx q[1];
rz(-2.8376104) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39897746) q[0];
sx q[0];
rz(-1.0053867) q[0];
sx q[0];
rz(-0.032755927) q[0];
rz(-pi) q[1];
rz(-1.4690222) q[2];
sx q[2];
rz(-1.8074805) q[2];
sx q[2];
rz(-0.19381154) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.099977255) q[1];
sx q[1];
rz(-1.4681889) q[1];
sx q[1];
rz(-2.2907738) q[1];
rz(2.899053) q[3];
sx q[3];
rz(-0.69787301) q[3];
sx q[3];
rz(-1.0353116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1659282) q[2];
sx q[2];
rz(-2.1037481) q[2];
sx q[2];
rz(0.69592875) q[2];
rz(1.1337229) q[3];
sx q[3];
rz(-1.1997831) q[3];
sx q[3];
rz(-2.5886562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7021779) q[0];
sx q[0];
rz(-1.004383) q[0];
sx q[0];
rz(-2.0942005) q[0];
rz(-1.9678496) q[1];
sx q[1];
rz(-2.4672697) q[1];
sx q[1];
rz(-1.6204999) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3701909) q[0];
sx q[0];
rz(-0.51734296) q[0];
sx q[0];
rz(-2.8964554) q[0];
rz(-0.32033605) q[2];
sx q[2];
rz(-0.89010677) q[2];
sx q[2];
rz(-2.5706511) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4075267) q[1];
sx q[1];
rz(-2.1868949) q[1];
sx q[1];
rz(-1.0259088) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71468227) q[3];
sx q[3];
rz(-0.41926786) q[3];
sx q[3];
rz(1.661232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.095470458) q[2];
sx q[2];
rz(-2.0746524) q[2];
sx q[2];
rz(2.6037237) q[2];
rz(1.6543903) q[3];
sx q[3];
rz(-2.7345149) q[3];
sx q[3];
rz(-2.4552104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9293514) q[0];
sx q[0];
rz(-0.23433267) q[0];
sx q[0];
rz(0.60337639) q[0];
rz(0.76250184) q[1];
sx q[1];
rz(-1.9489894) q[1];
sx q[1];
rz(2.4286483) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6978074) q[0];
sx q[0];
rz(-1.968211) q[0];
sx q[0];
rz(0.071594588) q[0];
x q[1];
rz(-1.845417) q[2];
sx q[2];
rz(-2.6634187) q[2];
sx q[2];
rz(2.630065) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7984287) q[1];
sx q[1];
rz(-2.6021932) q[1];
sx q[1];
rz(-1.1718962) q[1];
rz(2.5239046) q[3];
sx q[3];
rz(-2.0632944) q[3];
sx q[3];
rz(2.4503051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.9731628) q[2];
sx q[2];
rz(-0.91732401) q[2];
sx q[2];
rz(-1.470835) q[2];
rz(0.51703185) q[3];
sx q[3];
rz(-1.0443338) q[3];
sx q[3];
rz(-1.8487336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1626728) q[0];
sx q[0];
rz(-0.51386583) q[0];
sx q[0];
rz(-2.5901219) q[0];
rz(0.035471352) q[1];
sx q[1];
rz(-1.1330117) q[1];
sx q[1];
rz(1.6112526) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96824232) q[0];
sx q[0];
rz(-0.98778546) q[0];
sx q[0];
rz(0.75847404) q[0];
rz(0.82357652) q[2];
sx q[2];
rz(-1.608143) q[2];
sx q[2];
rz(2.3970791) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1332273) q[1];
sx q[1];
rz(-2.3043345) q[1];
sx q[1];
rz(-2.0627229) q[1];
rz(-pi) q[2];
rz(3.0569791) q[3];
sx q[3];
rz(-0.91949465) q[3];
sx q[3];
rz(-0.29825975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0772721) q[2];
sx q[2];
rz(-0.52017009) q[2];
sx q[2];
rz(1.2426097) q[2];
rz(0.51314917) q[3];
sx q[3];
rz(-2.7276701) q[3];
sx q[3];
rz(-1.7665524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.649491) q[0];
sx q[0];
rz(-0.92925564) q[0];
sx q[0];
rz(0.11665601) q[0];
rz(-0.77313441) q[1];
sx q[1];
rz(-1.1583068) q[1];
sx q[1];
rz(1.6366417) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3706544) q[0];
sx q[0];
rz(-2.2794834) q[0];
sx q[0];
rz(-0.38972008) q[0];
rz(-pi) q[1];
rz(-1.4261521) q[2];
sx q[2];
rz(-2.7491315) q[2];
sx q[2];
rz(1.5953428) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6500068) q[1];
sx q[1];
rz(-1.2131696) q[1];
sx q[1];
rz(-2.3050344) q[1];
rz(2.4890763) q[3];
sx q[3];
rz(-1.8015141) q[3];
sx q[3];
rz(2.0215061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2152805) q[2];
sx q[2];
rz(-0.24682385) q[2];
sx q[2];
rz(-0.8872633) q[2];
rz(-0.67982802) q[3];
sx q[3];
rz(-2.4726548) q[3];
sx q[3];
rz(0.63290709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55650869) q[0];
sx q[0];
rz(-1.8483138) q[0];
sx q[0];
rz(2.4853117) q[0];
rz(-2.9346924) q[1];
sx q[1];
rz(-1.2196536) q[1];
sx q[1];
rz(-1.9237178) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6164056) q[0];
sx q[0];
rz(-1.6058679) q[0];
sx q[0];
rz(0.082790815) q[0];
rz(2.0711259) q[2];
sx q[2];
rz(-2.4917951) q[2];
sx q[2];
rz(-1.4184679) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.52384842) q[1];
sx q[1];
rz(-0.19568014) q[1];
sx q[1];
rz(-0.62472549) q[1];
x q[2];
rz(2.7662606) q[3];
sx q[3];
rz(-2.4960244) q[3];
sx q[3];
rz(2.4433492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3394341) q[2];
sx q[2];
rz(-1.755244) q[2];
sx q[2];
rz(2.459724) q[2];
rz(-1.1784461) q[3];
sx q[3];
rz(-0.75740564) q[3];
sx q[3];
rz(1.2371548) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39663974) q[0];
sx q[0];
rz(-1.5322026) q[0];
sx q[0];
rz(2.9314281) q[0];
rz(2.919803) q[1];
sx q[1];
rz(-2.2237325) q[1];
sx q[1];
rz(0.04235696) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51084679) q[0];
sx q[0];
rz(-1.5486985) q[0];
sx q[0];
rz(2.9539897) q[0];
x q[1];
rz(0.63318166) q[2];
sx q[2];
rz(-2.0821619) q[2];
sx q[2];
rz(2.2610469) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.48407979) q[1];
sx q[1];
rz(-1.8575605) q[1];
sx q[1];
rz(0.57876719) q[1];
rz(-pi) q[2];
rz(0.70019763) q[3];
sx q[3];
rz(-2.4841016) q[3];
sx q[3];
rz(1.0669062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.28045851) q[2];
sx q[2];
rz(-2.0529842) q[2];
sx q[2];
rz(-2.4750278) q[2];
rz(2.376453) q[3];
sx q[3];
rz(-2.9355526) q[3];
sx q[3];
rz(-1.7407181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53369451) q[0];
sx q[0];
rz(-2.7602637) q[0];
sx q[0];
rz(-0.68699849) q[0];
rz(-1.2835361) q[1];
sx q[1];
rz(-1.1524408) q[1];
sx q[1];
rz(2.5680465) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4939369) q[0];
sx q[0];
rz(-1.0595317) q[0];
sx q[0];
rz(0.6880811) q[0];
x q[1];
rz(-0.23172371) q[2];
sx q[2];
rz(-2.0123693) q[2];
sx q[2];
rz(1.148846) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0189218) q[1];
sx q[1];
rz(-1.9098567) q[1];
sx q[1];
rz(2.1200722) q[1];
x q[2];
rz(2.403026) q[3];
sx q[3];
rz(-2.1436678) q[3];
sx q[3];
rz(2.3923739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0987229) q[2];
sx q[2];
rz(-1.9645773) q[2];
sx q[2];
rz(3.0549808) q[2];
rz(3.0232271) q[3];
sx q[3];
rz(-1.5990853) q[3];
sx q[3];
rz(-0.043244403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8260228) q[0];
sx q[0];
rz(-1.0181027) q[0];
sx q[0];
rz(-2.3636567) q[0];
rz(-2.8553873) q[1];
sx q[1];
rz(-0.83120167) q[1];
sx q[1];
rz(1.3298159) q[1];
rz(-1.329245) q[2];
sx q[2];
rz(-1.8941244) q[2];
sx q[2];
rz(2.6461442) q[2];
rz(0.56747464) q[3];
sx q[3];
rz(-2.6226061) q[3];
sx q[3];
rz(-2.0593986) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
