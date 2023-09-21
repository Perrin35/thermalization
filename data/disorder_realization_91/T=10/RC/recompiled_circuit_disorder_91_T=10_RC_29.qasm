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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5509697) q[0];
sx q[0];
rz(-1.5184214) q[0];
sx q[0];
rz(2.4916324) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6257329) q[2];
sx q[2];
rz(-1.5896279) q[2];
sx q[2];
rz(0.042382391) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6273856) q[1];
sx q[1];
rz(-1.5878907) q[1];
sx q[1];
rz(1.8219201) q[1];
rz(-pi) q[2];
rz(-2.8475548) q[3];
sx q[3];
rz(-0.66411823) q[3];
sx q[3];
rz(-0.77120632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3657637) q[2];
sx q[2];
rz(-2.3684431) q[2];
sx q[2];
rz(-1.3295056) q[2];
rz(1.4154411) q[3];
sx q[3];
rz(-1.5751782) q[3];
sx q[3];
rz(0.20234385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99130327) q[0];
sx q[0];
rz(-0.54420272) q[0];
sx q[0];
rz(1.4527028) q[0];
rz(2.9406722) q[1];
sx q[1];
rz(-2.0566437) q[1];
sx q[1];
rz(-0.30028775) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9584504) q[0];
sx q[0];
rz(-2.169325) q[0];
sx q[0];
rz(-1.7167702) q[0];
rz(-pi) q[1];
rz(1.8882206) q[2];
sx q[2];
rz(-1.023264) q[2];
sx q[2];
rz(-0.65371338) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.049388) q[1];
sx q[1];
rz(-0.37282473) q[1];
sx q[1];
rz(-0.91918175) q[1];
rz(1.3817915) q[3];
sx q[3];
rz(-1.4100037) q[3];
sx q[3];
rz(-2.8652428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.369027) q[2];
sx q[2];
rz(-0.3521266) q[2];
sx q[2];
rz(1.0380113) q[2];
rz(0.84233061) q[3];
sx q[3];
rz(-1.4146283) q[3];
sx q[3];
rz(2.7712908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5623986) q[0];
sx q[0];
rz(-2.6694522) q[0];
sx q[0];
rz(2.753479) q[0];
rz(0.072470486) q[1];
sx q[1];
rz(-1.4265172) q[1];
sx q[1];
rz(-2.8252576) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4916723) q[0];
sx q[0];
rz(-0.24501093) q[0];
sx q[0];
rz(1.5632167) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5691368) q[2];
sx q[2];
rz(-1.2081895) q[2];
sx q[2];
rz(0.87871694) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9005147) q[1];
sx q[1];
rz(-1.9652275) q[1];
sx q[1];
rz(2.6532252) q[1];
rz(-pi) q[2];
rz(-0.29052492) q[3];
sx q[3];
rz(-1.2447262) q[3];
sx q[3];
rz(-0.80929221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1091653) q[2];
sx q[2];
rz(-2.9563603) q[2];
sx q[2];
rz(-0.16080984) q[2];
rz(3.0155449) q[3];
sx q[3];
rz(-1.3879317) q[3];
sx q[3];
rz(1.0236615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9008824) q[0];
sx q[0];
rz(-2.4375589) q[0];
sx q[0];
rz(2.3098992) q[0];
rz(-1.7968934) q[1];
sx q[1];
rz(-0.99935499) q[1];
sx q[1];
rz(1.6279189) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3095703) q[0];
sx q[0];
rz(-2.4674468) q[0];
sx q[0];
rz(-0.073622965) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.081110031) q[2];
sx q[2];
rz(-1.0081648) q[2];
sx q[2];
rz(-2.5155591) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0107683) q[1];
sx q[1];
rz(-1.1225892) q[1];
sx q[1];
rz(0.018647714) q[1];
rz(-2.4098445) q[3];
sx q[3];
rz(-1.2536612) q[3];
sx q[3];
rz(1.5750615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4251129) q[2];
sx q[2];
rz(-2.1322865) q[2];
sx q[2];
rz(1.224219) q[2];
rz(2.7246357) q[3];
sx q[3];
rz(-1.0427534) q[3];
sx q[3];
rz(2.6127889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083374627) q[0];
sx q[0];
rz(-2.871802) q[0];
sx q[0];
rz(1.8378687) q[0];
rz(-0.45571348) q[1];
sx q[1];
rz(-0.2625176) q[1];
sx q[1];
rz(-0.051503332) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18729678) q[0];
sx q[0];
rz(-1.4547537) q[0];
sx q[0];
rz(0.10033484) q[0];
rz(1.2850045) q[2];
sx q[2];
rz(-1.7190061) q[2];
sx q[2];
rz(-2.1697901) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0806549) q[1];
sx q[1];
rz(-2.6662711) q[1];
sx q[1];
rz(1.2259543) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6938985) q[3];
sx q[3];
rz(-0.8408635) q[3];
sx q[3];
rz(-0.52809944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6872528) q[2];
sx q[2];
rz(-1.3663102) q[2];
sx q[2];
rz(-2.5435737) q[2];
rz(0.55001843) q[3];
sx q[3];
rz(-0.23967448) q[3];
sx q[3];
rz(-2.0050744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0710058) q[0];
sx q[0];
rz(-3.0650009) q[0];
sx q[0];
rz(0.20198527) q[0];
rz(0.96549353) q[1];
sx q[1];
rz(-1.0243203) q[1];
sx q[1];
rz(-3.0583256) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71654728) q[0];
sx q[0];
rz(-1.3971359) q[0];
sx q[0];
rz(0.19787775) q[0];
rz(-pi) q[1];
rz(1.303057) q[2];
sx q[2];
rz(-1.8362852) q[2];
sx q[2];
rz(-2.7051085) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6492918) q[1];
sx q[1];
rz(-1.9363083) q[1];
sx q[1];
rz(-3.027012) q[1];
rz(0.55732754) q[3];
sx q[3];
rz(-1.8263655) q[3];
sx q[3];
rz(-2.717358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.10963708) q[2];
sx q[2];
rz(-1.457931) q[2];
sx q[2];
rz(1.3827682) q[2];
rz(-2.8921195) q[3];
sx q[3];
rz(-1.8981257) q[3];
sx q[3];
rz(-1.4613387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7145342) q[0];
sx q[0];
rz(-2.336851) q[0];
sx q[0];
rz(-0.89685857) q[0];
rz(-2.8915021) q[1];
sx q[1];
rz(-0.18083328) q[1];
sx q[1];
rz(2.0239963) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9761336) q[0];
sx q[0];
rz(-1.4384067) q[0];
sx q[0];
rz(2.4112941) q[0];
rz(-pi) q[1];
rz(-1.8115225) q[2];
sx q[2];
rz(-2.4154818) q[2];
sx q[2];
rz(-1.7118529) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7270131) q[1];
sx q[1];
rz(-2.1257183) q[1];
sx q[1];
rz(1.8068061) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.696633) q[3];
sx q[3];
rz(-2.5611454) q[3];
sx q[3];
rz(2.2821033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.08427944) q[2];
sx q[2];
rz(-1.3998312) q[2];
sx q[2];
rz(0.21952595) q[2];
rz(-0.38671842) q[3];
sx q[3];
rz(-2.3733449) q[3];
sx q[3];
rz(-0.99532551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052208386) q[0];
sx q[0];
rz(-1.5445671) q[0];
sx q[0];
rz(-0.21729939) q[0];
rz(-0.030933881) q[1];
sx q[1];
rz(-0.63806454) q[1];
sx q[1];
rz(2.4826179) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9650759) q[0];
sx q[0];
rz(-0.79916164) q[0];
sx q[0];
rz(1.6060711) q[0];
rz(-pi) q[1];
rz(-0.73924139) q[2];
sx q[2];
rz(-0.87202245) q[2];
sx q[2];
rz(-0.07585635) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36775667) q[1];
sx q[1];
rz(-0.78513297) q[1];
sx q[1];
rz(0.30898856) q[1];
rz(-pi) q[2];
rz(0.92774763) q[3];
sx q[3];
rz(-1.5156931) q[3];
sx q[3];
rz(2.7411214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4387681) q[2];
sx q[2];
rz(-1.7610901) q[2];
sx q[2];
rz(2.8213815) q[2];
rz(-1.5252339) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(-1.2493791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76072389) q[0];
sx q[0];
rz(-2.9601233) q[0];
sx q[0];
rz(2.4587801) q[0];
rz(2.071351) q[1];
sx q[1];
rz(-1.2425334) q[1];
sx q[1];
rz(1.210093) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34459201) q[0];
sx q[0];
rz(-1.7779777) q[0];
sx q[0];
rz(1.209757) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8898583) q[2];
sx q[2];
rz(-0.98099698) q[2];
sx q[2];
rz(-2.7272607) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.71483892) q[1];
sx q[1];
rz(-2.4086082) q[1];
sx q[1];
rz(1.1232125) q[1];
x q[2];
rz(0.61981598) q[3];
sx q[3];
rz(-1.2634988) q[3];
sx q[3];
rz(0.67051552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.575763) q[2];
sx q[2];
rz(-2.1255707) q[2];
sx q[2];
rz(0.56662095) q[2];
rz(2.2120655) q[3];
sx q[3];
rz(-0.98171392) q[3];
sx q[3];
rz(1.367759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.728445) q[0];
sx q[0];
rz(-2.8840273) q[0];
sx q[0];
rz(1.7096827) q[0];
rz(0.52629772) q[1];
sx q[1];
rz(-0.52733517) q[1];
sx q[1];
rz(2.3419103) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47213263) q[0];
sx q[0];
rz(-1.2233943) q[0];
sx q[0];
rz(-3.0513289) q[0];
rz(-pi) q[1];
rz(0.16074796) q[2];
sx q[2];
rz(-1.4904009) q[2];
sx q[2];
rz(2.423167) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4754776) q[1];
sx q[1];
rz(-1.0628504) q[1];
sx q[1];
rz(0.39132262) q[1];
rz(-pi) q[2];
rz(-0.14567104) q[3];
sx q[3];
rz(-1.4958069) q[3];
sx q[3];
rz(-2.2742227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8563103) q[2];
sx q[2];
rz(-0.44390634) q[2];
sx q[2];
rz(-0.69520673) q[2];
rz(-0.55784145) q[3];
sx q[3];
rz(-1.3643684) q[3];
sx q[3];
rz(0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.0556864) q[0];
sx q[0];
rz(-1.9491371) q[0];
sx q[0];
rz(0.51167713) q[0];
rz(1.2790537) q[1];
sx q[1];
rz(-0.74395724) q[1];
sx q[1];
rz(-0.67768135) q[1];
rz(3.1391115) q[2];
sx q[2];
rz(-2.8056792) q[2];
sx q[2];
rz(0.16028595) q[2];
rz(1.4064895) q[3];
sx q[3];
rz(-0.77948979) q[3];
sx q[3];
rz(-0.43850552) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];