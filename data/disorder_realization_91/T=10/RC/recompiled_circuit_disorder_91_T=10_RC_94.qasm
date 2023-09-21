OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.70541731) q[0];
sx q[0];
rz(-2.5751312) q[0];
sx q[0];
rz(-0.17106549) q[0];
rz(-1.6879727) q[1];
sx q[1];
rz(-2.7984518) q[1];
sx q[1];
rz(-1.8106102) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5509697) q[0];
sx q[0];
rz(-1.5184214) q[0];
sx q[0];
rz(2.4916324) q[0];
x q[1];
rz(1.5491484) q[2];
sx q[2];
rz(-2.0865555) q[2];
sx q[2];
rz(-1.6025008) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.060974412) q[1];
sx q[1];
rz(-1.3197101) q[1];
sx q[1];
rz(3.1239448) q[1];
rz(0.29403789) q[3];
sx q[3];
rz(-2.4774744) q[3];
sx q[3];
rz(0.77120632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.77582899) q[2];
sx q[2];
rz(-0.77314955) q[2];
sx q[2];
rz(1.3295056) q[2];
rz(-1.4154411) q[3];
sx q[3];
rz(-1.5664145) q[3];
sx q[3];
rz(0.20234385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1502894) q[0];
sx q[0];
rz(-0.54420272) q[0];
sx q[0];
rz(-1.6888899) q[0];
rz(2.9406722) q[1];
sx q[1];
rz(-2.0566437) q[1];
sx q[1];
rz(-0.30028775) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0694885) q[0];
sx q[0];
rz(-2.527643) q[0];
sx q[0];
rz(-2.9314562) q[0];
rz(-0.57057256) q[2];
sx q[2];
rz(-1.3010446) q[2];
sx q[2];
rz(-2.3938993) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2806432) q[1];
sx q[1];
rz(-1.7935392) q[1];
sx q[1];
rz(-1.8722948) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85864752) q[3];
sx q[3];
rz(-0.24752366) q[3];
sx q[3];
rz(-1.1502707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.369027) q[2];
sx q[2];
rz(-2.7894661) q[2];
sx q[2];
rz(-1.0380113) q[2];
rz(2.299262) q[3];
sx q[3];
rz(-1.4146283) q[3];
sx q[3];
rz(0.37030181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(1.5623986) q[0];
sx q[0];
rz(-2.6694522) q[0];
sx q[0];
rz(0.38811362) q[0];
rz(3.0691222) q[1];
sx q[1];
rz(-1.7150755) q[1];
sx q[1];
rz(0.31633502) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2280699) q[0];
sx q[0];
rz(-1.5689578) q[0];
sx q[0];
rz(-1.3257922) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5306273) q[2];
sx q[2];
rz(-0.66662153) q[2];
sx q[2];
rz(1.9463469) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2410779) q[1];
sx q[1];
rz(-1.1763651) q[1];
sx q[1];
rz(-2.6532252) q[1];
rz(0.86795904) q[3];
sx q[3];
rz(-0.43324019) q[3];
sx q[3];
rz(0.058335282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0324273) q[2];
sx q[2];
rz(-2.9563603) q[2];
sx q[2];
rz(0.16080984) q[2];
rz(0.12604776) q[3];
sx q[3];
rz(-1.3879317) q[3];
sx q[3];
rz(2.1179312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2407103) q[0];
sx q[0];
rz(-2.4375589) q[0];
sx q[0];
rz(-2.3098992) q[0];
rz(-1.7968934) q[1];
sx q[1];
rz(-0.99935499) q[1];
sx q[1];
rz(-1.5136738) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31878372) q[0];
sx q[0];
rz(-1.6167287) q[0];
sx q[0];
rz(0.6728234) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6985745) q[2];
sx q[2];
rz(-2.5737692) q[2];
sx q[2];
rz(0.77726269) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0537783) q[1];
sx q[1];
rz(-2.6930241) q[1];
sx q[1];
rz(-1.5320369) q[1];
rz(-pi) q[2];
rz(1.1553331) q[3];
sx q[3];
rz(-2.2586125) q[3];
sx q[3];
rz(-2.8727939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4251129) q[2];
sx q[2];
rz(-1.0093062) q[2];
sx q[2];
rz(-1.9173737) q[2];
rz(0.41695693) q[3];
sx q[3];
rz(-2.0988393) q[3];
sx q[3];
rz(-0.52880374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.083374627) q[0];
sx q[0];
rz(-2.871802) q[0];
sx q[0];
rz(-1.303724) q[0];
rz(0.45571348) q[1];
sx q[1];
rz(-2.8790751) q[1];
sx q[1];
rz(-0.051503332) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6131825) q[0];
sx q[0];
rz(-2.9883356) q[0];
sx q[0];
rz(2.2806703) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2850045) q[2];
sx q[2];
rz(-1.7190061) q[2];
sx q[2];
rz(2.1697901) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4645849) q[1];
sx q[1];
rz(-2.0160463) q[1];
sx q[1];
rz(-2.9693309) q[1];
rz(-pi) q[2];
rz(3.0052161) q[3];
sx q[3];
rz(-0.73835056) q[3];
sx q[3];
rz(-2.7969558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6872528) q[2];
sx q[2];
rz(-1.7752825) q[2];
sx q[2];
rz(-2.5435737) q[2];
rz(2.5915742) q[3];
sx q[3];
rz(-2.9019182) q[3];
sx q[3];
rz(1.1365183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0705868) q[0];
sx q[0];
rz(-3.0650009) q[0];
sx q[0];
rz(-0.20198527) q[0];
rz(-0.96549353) q[1];
sx q[1];
rz(-1.0243203) q[1];
sx q[1];
rz(3.0583256) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9986388) q[0];
sx q[0];
rz(-2.8790701) q[0];
sx q[0];
rz(-0.72857626) q[0];
rz(0.27481885) q[2];
sx q[2];
rz(-1.3126557) q[2];
sx q[2];
rz(-1.9354265) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0373842) q[1];
sx q[1];
rz(-1.4638149) q[1];
sx q[1];
rz(1.9385098) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5842651) q[3];
sx q[3];
rz(-1.3152272) q[3];
sx q[3];
rz(-0.42423466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0319556) q[2];
sx q[2];
rz(-1.457931) q[2];
sx q[2];
rz(1.7588245) q[2];
rz(2.8921195) q[3];
sx q[3];
rz(-1.8981257) q[3];
sx q[3];
rz(-1.680254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7145342) q[0];
sx q[0];
rz(-0.80474168) q[0];
sx q[0];
rz(0.89685857) q[0];
rz(-0.25009051) q[1];
sx q[1];
rz(-0.18083328) q[1];
sx q[1];
rz(1.1175964) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5516978) q[0];
sx q[0];
rz(-2.4015744) q[0];
sx q[0];
rz(0.19703534) q[0];
rz(-pi) q[1];
x q[1];
rz(0.859185) q[2];
sx q[2];
rz(-1.4118328) q[2];
sx q[2];
rz(0.32260103) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0301789) q[1];
sx q[1];
rz(-1.3707268) q[1];
sx q[1];
rz(2.5740037) q[1];
rz(-pi) q[2];
rz(-2.696633) q[3];
sx q[3];
rz(-0.58044725) q[3];
sx q[3];
rz(-2.2821033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.08427944) q[2];
sx q[2];
rz(-1.3998312) q[2];
sx q[2];
rz(-0.21952595) q[2];
rz(0.38671842) q[3];
sx q[3];
rz(-2.3733449) q[3];
sx q[3];
rz(-2.1462671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(0.052208386) q[0];
sx q[0];
rz(-1.5445671) q[0];
sx q[0];
rz(-2.9242933) q[0];
rz(-3.1106588) q[1];
sx q[1];
rz(-2.5035281) q[1];
sx q[1];
rz(-0.6589748) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17651672) q[0];
sx q[0];
rz(-0.79916164) q[0];
sx q[0];
rz(-1.6060711) q[0];
rz(-pi) q[1];
rz(0.73924139) q[2];
sx q[2];
rz(-0.87202245) q[2];
sx q[2];
rz(-3.0657363) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.79170376) q[1];
sx q[1];
rz(-0.83194299) q[1];
sx q[1];
rz(-1.2757343) q[1];
x q[2];
rz(0.92774763) q[3];
sx q[3];
rz(-1.5156931) q[3];
sx q[3];
rz(2.7411214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4387681) q[2];
sx q[2];
rz(-1.3805026) q[2];
sx q[2];
rz(0.32021114) q[2];
rz(1.6163588) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(1.8922136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3808688) q[0];
sx q[0];
rz(-0.18146935) q[0];
sx q[0];
rz(0.6828126) q[0];
rz(1.0702417) q[1];
sx q[1];
rz(-1.2425334) q[1];
sx q[1];
rz(-1.210093) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34459201) q[0];
sx q[0];
rz(-1.3636149) q[0];
sx q[0];
rz(-1.209757) q[0];
x q[1];
rz(-1.2145043) q[2];
sx q[2];
rz(-0.63535832) q[2];
sx q[2];
rz(-2.2941342) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.51325071) q[1];
sx q[1];
rz(-1.2770137) q[1];
sx q[1];
rz(-2.252584) q[1];
rz(-pi) q[2];
rz(-2.5217767) q[3];
sx q[3];
rz(-1.2634988) q[3];
sx q[3];
rz(-2.4710771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5658297) q[2];
sx q[2];
rz(-1.016022) q[2];
sx q[2];
rz(-2.5749717) q[2];
rz(-0.9295272) q[3];
sx q[3];
rz(-2.1598787) q[3];
sx q[3];
rz(-1.367759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.728445) q[0];
sx q[0];
rz(-2.8840273) q[0];
sx q[0];
rz(1.7096827) q[0];
rz(2.6152949) q[1];
sx q[1];
rz(-2.6142575) q[1];
sx q[1];
rz(-0.79968232) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47213263) q[0];
sx q[0];
rz(-1.2233943) q[0];
sx q[0];
rz(-3.0513289) q[0];
rz(-pi) q[1];
rz(1.6522371) q[2];
sx q[2];
rz(-1.7310206) q[2];
sx q[2];
rz(-2.3022431) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.77196808) q[1];
sx q[1];
rz(-2.5110285) q[1];
sx q[1];
rz(-0.97009138) q[1];
rz(-pi) q[2];
rz(-0.14567104) q[3];
sx q[3];
rz(-1.6457857) q[3];
sx q[3];
rz(-0.86736995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8563103) q[2];
sx q[2];
rz(-2.6976863) q[2];
sx q[2];
rz(2.4463859) q[2];
rz(-2.5837512) q[3];
sx q[3];
rz(-1.7772243) q[3];
sx q[3];
rz(0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0556864) q[0];
sx q[0];
rz(-1.9491371) q[0];
sx q[0];
rz(0.51167713) q[0];
rz(1.862539) q[1];
sx q[1];
rz(-2.3976354) q[1];
sx q[1];
rz(2.4639113) q[1];
rz(2.8056801) q[2];
sx q[2];
rz(-1.5716142) q[2];
sx q[2];
rz(-1.4128528) q[2];
rz(0.79808509) q[3];
sx q[3];
rz(-1.4555664) q[3];
sx q[3];
rz(1.2496787) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
