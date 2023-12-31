OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4296071) q[0];
sx q[0];
rz(-0.40661943) q[0];
sx q[0];
rz(-2.8924195) q[0];
rz(3.0781526) q[1];
sx q[1];
rz(-0.97172207) q[1];
sx q[1];
rz(2.5914153) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8367856) q[0];
sx q[0];
rz(-1.7903622) q[0];
sx q[0];
rz(0.027713393) q[0];
rz(-2.1044188) q[2];
sx q[2];
rz(-1.954477) q[2];
sx q[2];
rz(-1.2984315) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.3783592) q[1];
sx q[1];
rz(-2.2765056) q[1];
sx q[1];
rz(2.6325429) q[1];
rz(2.1816741) q[3];
sx q[3];
rz(-1.6392518) q[3];
sx q[3];
rz(-1.6144891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.41574079) q[2];
sx q[2];
rz(-0.44885138) q[2];
sx q[2];
rz(-0.63981167) q[2];
rz(0.85302991) q[3];
sx q[3];
rz(-2.5363688) q[3];
sx q[3];
rz(-0.38133347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8752276) q[0];
sx q[0];
rz(-1.9748283) q[0];
sx q[0];
rz(-0.27045989) q[0];
rz(-2.4282783) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(-1.5126022) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0727901) q[0];
sx q[0];
rz(-1.8913942) q[0];
sx q[0];
rz(-2.1945303) q[0];
rz(3.1299635) q[2];
sx q[2];
rz(-2.8577869) q[2];
sx q[2];
rz(-2.0276045) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4334129) q[1];
sx q[1];
rz(-0.69060329) q[1];
sx q[1];
rz(1.3604926) q[1];
rz(-pi) q[2];
rz(-0.87725957) q[3];
sx q[3];
rz(-2.7765397) q[3];
sx q[3];
rz(-3.0960494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2018532) q[2];
sx q[2];
rz(-0.24213232) q[2];
sx q[2];
rz(-0.80336037) q[2];
rz(-1.057829) q[3];
sx q[3];
rz(-1.4927031) q[3];
sx q[3];
rz(0.025618205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8818883) q[0];
sx q[0];
rz(-2.0443679) q[0];
sx q[0];
rz(-0.29552466) q[0];
rz(0.23513901) q[1];
sx q[1];
rz(-1.7087015) q[1];
sx q[1];
rz(0.74584109) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5416109) q[0];
sx q[0];
rz(-2.9020502) q[0];
sx q[0];
rz(1.9639067) q[0];
x q[1];
rz(-2.8638641) q[2];
sx q[2];
rz(-0.54364294) q[2];
sx q[2];
rz(0.03253983) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5513788) q[1];
sx q[1];
rz(-2.5775902) q[1];
sx q[1];
rz(2.6022634) q[1];
x q[2];
rz(2.6607473) q[3];
sx q[3];
rz(-0.50243176) q[3];
sx q[3];
rz(-1.06711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.088034257) q[2];
sx q[2];
rz(-0.025667889) q[2];
sx q[2];
rz(-2.4528465) q[2];
rz(-3.0912494) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(-1.6200199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21928366) q[0];
sx q[0];
rz(-2.121121) q[0];
sx q[0];
rz(0.10198378) q[0];
rz(0.12022262) q[1];
sx q[1];
rz(-2.6133803) q[1];
sx q[1];
rz(0.27329683) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2647588) q[0];
sx q[0];
rz(-2.8059373) q[0];
sx q[0];
rz(-1.8434974) q[0];
rz(-0.64812135) q[2];
sx q[2];
rz(-1.1037877) q[2];
sx q[2];
rz(-1.8304706) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3015775) q[1];
sx q[1];
rz(-1.8043648) q[1];
sx q[1];
rz(2.9905) q[1];
x q[2];
rz(-0.23394211) q[3];
sx q[3];
rz(-2.4787239) q[3];
sx q[3];
rz(2.5556504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3499202) q[2];
sx q[2];
rz(-0.73799729) q[2];
sx q[2];
rz(0.50393528) q[2];
rz(0.079581633) q[3];
sx q[3];
rz(-1.1527529) q[3];
sx q[3];
rz(2.7664405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85916096) q[0];
sx q[0];
rz(-2.0895884) q[0];
sx q[0];
rz(3.058847) q[0];
rz(0.67963183) q[1];
sx q[1];
rz(-1.6487164) q[1];
sx q[1];
rz(-2.1544429) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3521096) q[0];
sx q[0];
rz(-1.9728567) q[0];
sx q[0];
rz(-0.91870086) q[0];
x q[1];
rz(0.64381386) q[2];
sx q[2];
rz(-1.5632731) q[2];
sx q[2];
rz(1.3155589) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.134569) q[1];
sx q[1];
rz(-0.94788523) q[1];
sx q[1];
rz(1.9593777) q[1];
rz(-pi) q[2];
rz(-1.9636642) q[3];
sx q[3];
rz(-1.53189) q[3];
sx q[3];
rz(-1.2272569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8905028) q[2];
sx q[2];
rz(-2.3513998) q[2];
sx q[2];
rz(-2.9303072) q[2];
rz(-2.7206897) q[3];
sx q[3];
rz(-0.55287164) q[3];
sx q[3];
rz(3.0781854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6565276) q[0];
sx q[0];
rz(-1.2542897) q[0];
sx q[0];
rz(-0.791839) q[0];
rz(-0.99545288) q[1];
sx q[1];
rz(-0.96770006) q[1];
sx q[1];
rz(1.8621559) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56023843) q[0];
sx q[0];
rz(-1.2890153) q[0];
sx q[0];
rz(-0.015334107) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60792578) q[2];
sx q[2];
rz(-1.8106917) q[2];
sx q[2];
rz(1.243967) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2710323) q[1];
sx q[1];
rz(-0.522627) q[1];
sx q[1];
rz(1.51103) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2313314) q[3];
sx q[3];
rz(-2.913774) q[3];
sx q[3];
rz(-1.4753301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.25097686) q[2];
sx q[2];
rz(-1.0303409) q[2];
sx q[2];
rz(2.3708169) q[2];
rz(-1.4701014) q[3];
sx q[3];
rz(-2.7272868) q[3];
sx q[3];
rz(2.9582086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.32271785) q[0];
sx q[0];
rz(-0.6495496) q[0];
sx q[0];
rz(3.0859257) q[0];
rz(0.21559134) q[1];
sx q[1];
rz(-0.76342738) q[1];
sx q[1];
rz(-0.0035704426) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5033747) q[0];
sx q[0];
rz(-0.57045454) q[0];
sx q[0];
rz(0.72871448) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3891719) q[2];
sx q[2];
rz(-1.2781029) q[2];
sx q[2];
rz(0.21586403) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.53686245) q[1];
sx q[1];
rz(-2.4559048) q[1];
sx q[1];
rz(2.0525949) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7690897) q[3];
sx q[3];
rz(-1.4962422) q[3];
sx q[3];
rz(-1.2411181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.59721649) q[2];
sx q[2];
rz(-0.95149136) q[2];
sx q[2];
rz(-0.92010951) q[2];
rz(0.19872935) q[3];
sx q[3];
rz(-1.2343497) q[3];
sx q[3];
rz(0.41771093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51628095) q[0];
sx q[0];
rz(-1.6315062) q[0];
sx q[0];
rz(1.4177119) q[0];
rz(-2.7334546) q[1];
sx q[1];
rz(-2.1022271) q[1];
sx q[1];
rz(-0.67869854) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0146101) q[0];
sx q[0];
rz(-1.8982366) q[0];
sx q[0];
rz(1.2922657) q[0];
rz(0.46531123) q[2];
sx q[2];
rz(-0.88854549) q[2];
sx q[2];
rz(-1.0516143) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.39290909) q[1];
sx q[1];
rz(-0.52782413) q[1];
sx q[1];
rz(-1.4455568) q[1];
rz(-2.802556) q[3];
sx q[3];
rz(-2.6936274) q[3];
sx q[3];
rz(0.26089222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9644908) q[2];
sx q[2];
rz(-1.0937546) q[2];
sx q[2];
rz(-0.91782451) q[2];
rz(-1.9994036) q[3];
sx q[3];
rz(-0.85421383) q[3];
sx q[3];
rz(-0.93723047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1552102) q[0];
sx q[0];
rz(-0.6739524) q[0];
sx q[0];
rz(-0.25979364) q[0];
rz(-0.70867509) q[1];
sx q[1];
rz(-2.8627113) q[1];
sx q[1];
rz(-0.52694595) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5848815) q[0];
sx q[0];
rz(-1.3538133) q[0];
sx q[0];
rz(-2.7816539) q[0];
rz(-0.91295816) q[2];
sx q[2];
rz(-1.3557938) q[2];
sx q[2];
rz(-1.022033) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.38339112) q[1];
sx q[1];
rz(-1.8806629) q[1];
sx q[1];
rz(1.8499225) q[1];
rz(-1.8606887) q[3];
sx q[3];
rz(-1.999951) q[3];
sx q[3];
rz(0.030062519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.32593411) q[2];
sx q[2];
rz(-0.84647536) q[2];
sx q[2];
rz(2.3507067) q[2];
rz(-0.38665006) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(-2.8044243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7336422) q[0];
sx q[0];
rz(-0.17245094) q[0];
sx q[0];
rz(-0.98544425) q[0];
rz(2.573029) q[1];
sx q[1];
rz(-1.111258) q[1];
sx q[1];
rz(-2.7808166) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4496778) q[0];
sx q[0];
rz(-1.5957818) q[0];
sx q[0];
rz(-1.7382394) q[0];
rz(-pi) q[1];
rz(-1.785727) q[2];
sx q[2];
rz(-1.4533918) q[2];
sx q[2];
rz(-1.6984775) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0204748) q[1];
sx q[1];
rz(-1.4282244) q[1];
sx q[1];
rz(-2.3974182) q[1];
x q[2];
rz(1.4233227) q[3];
sx q[3];
rz(-1.7383988) q[3];
sx q[3];
rz(2.9180805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0439904) q[2];
sx q[2];
rz(-2.4520935) q[2];
sx q[2];
rz(-2.1981751) q[2];
rz(-0.57389456) q[3];
sx q[3];
rz(-2.6608163) q[3];
sx q[3];
rz(-1.3963612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5744793) q[0];
sx q[0];
rz(-1.4470826) q[0];
sx q[0];
rz(-0.8599109) q[0];
rz(1.7815331) q[1];
sx q[1];
rz(-2.139745) q[1];
sx q[1];
rz(2.3812961) q[1];
rz(-3.0439539) q[2];
sx q[2];
rz(-1.7751481) q[2];
sx q[2];
rz(-2.912599) q[2];
rz(0.73137024) q[3];
sx q[3];
rz(-0.89805713) q[3];
sx q[3];
rz(-3.0039207) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
