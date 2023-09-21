OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1155137) q[0];
sx q[0];
rz(-1.4839412) q[0];
sx q[0];
rz(2.8154362) q[0];
rz(-1.1905319) q[1];
sx q[1];
rz(-1.3500554) q[1];
sx q[1];
rz(-1.5989369) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9261949) q[0];
sx q[0];
rz(-2.8370259) q[0];
sx q[0];
rz(-2.347441) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5390981) q[2];
sx q[2];
rz(-1.3598816) q[2];
sx q[2];
rz(-0.22533016) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.76347683) q[1];
sx q[1];
rz(-2.8045178) q[1];
sx q[1];
rz(2.377305) q[1];
rz(-1.5845675) q[3];
sx q[3];
rz(-0.78409401) q[3];
sx q[3];
rz(2.9053094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8618384) q[2];
sx q[2];
rz(-2.162231) q[2];
sx q[2];
rz(-2.2564783) q[2];
rz(2.4195813) q[3];
sx q[3];
rz(-1.4530028) q[3];
sx q[3];
rz(0.0074145934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9279813) q[0];
sx q[0];
rz(-0.95887029) q[0];
sx q[0];
rz(-1.0990748) q[0];
rz(0.66501578) q[1];
sx q[1];
rz(-1.4140833) q[1];
sx q[1];
rz(2.2639993) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7769593) q[0];
sx q[0];
rz(-1.3898802) q[0];
sx q[0];
rz(-0.69676708) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7769233) q[2];
sx q[2];
rz(-1.696535) q[2];
sx q[2];
rz(-2.113935) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6192012) q[1];
sx q[1];
rz(-0.11208216) q[1];
sx q[1];
rz(-1.2680608) q[1];
rz(-0.86032805) q[3];
sx q[3];
rz(-0.48947696) q[3];
sx q[3];
rz(1.5599172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.39891222) q[2];
sx q[2];
rz(-1.8360527) q[2];
sx q[2];
rz(-1.1304643) q[2];
rz(1.8418664) q[3];
sx q[3];
rz(-1.9176509) q[3];
sx q[3];
rz(-1.4484423) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.995342) q[0];
sx q[0];
rz(-1.2820219) q[0];
sx q[0];
rz(-2.8515942) q[0];
rz(-2.4747804) q[1];
sx q[1];
rz(-1.0338444) q[1];
sx q[1];
rz(0.07382948) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8683257) q[0];
sx q[0];
rz(-3.0330015) q[0];
sx q[0];
rz(-0.6032087) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29088144) q[2];
sx q[2];
rz(-1.8996432) q[2];
sx q[2];
rz(1.8287303) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.044588305) q[1];
sx q[1];
rz(-1.6267596) q[1];
sx q[1];
rz(1.7486649) q[1];
x q[2];
rz(-1.4521493) q[3];
sx q[3];
rz(-1.2216179) q[3];
sx q[3];
rz(0.94482869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4905711) q[2];
sx q[2];
rz(-1.1940424) q[2];
sx q[2];
rz(-2.4948965) q[2];
rz(-1.1086639) q[3];
sx q[3];
rz(-0.78287786) q[3];
sx q[3];
rz(-1.1289319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9639503) q[0];
sx q[0];
rz(-0.17340604) q[0];
sx q[0];
rz(1.9529163) q[0];
rz(-1.0186609) q[1];
sx q[1];
rz(-0.97266346) q[1];
sx q[1];
rz(1.7046938) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80285145) q[0];
sx q[0];
rz(-0.84206284) q[0];
sx q[0];
rz(-1.4472423) q[0];
rz(1.2595348) q[2];
sx q[2];
rz(-1.2584104) q[2];
sx q[2];
rz(-0.54158467) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.97836271) q[1];
sx q[1];
rz(-2.0320738) q[1];
sx q[1];
rz(0.39050885) q[1];
rz(-pi) q[2];
rz(1.7014245) q[3];
sx q[3];
rz(-1.4454953) q[3];
sx q[3];
rz(0.13392042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.556095) q[2];
sx q[2];
rz(-1.4336339) q[2];
sx q[2];
rz(-2.0193224) q[2];
rz(1.026011) q[3];
sx q[3];
rz(-2.3882073) q[3];
sx q[3];
rz(-0.99075738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4595903) q[0];
sx q[0];
rz(-0.90072173) q[0];
sx q[0];
rz(0.37297747) q[0];
rz(-2.9176118) q[1];
sx q[1];
rz(-1.9517027) q[1];
sx q[1];
rz(-1.3164828) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7861159) q[0];
sx q[0];
rz(-1.9516203) q[0];
sx q[0];
rz(-2.064408) q[0];
rz(-pi) q[1];
rz(-2.3773642) q[2];
sx q[2];
rz(-2.7856305) q[2];
sx q[2];
rz(2.228235) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.59606325) q[1];
sx q[1];
rz(-1.5585594) q[1];
sx q[1];
rz(1.5188602) q[1];
x q[2];
rz(-2.5943807) q[3];
sx q[3];
rz(-2.0378761) q[3];
sx q[3];
rz(-0.78391778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0405154) q[2];
sx q[2];
rz(-0.83156362) q[2];
sx q[2];
rz(0.80580795) q[2];
rz(-0.53330437) q[3];
sx q[3];
rz(-1.1321944) q[3];
sx q[3];
rz(-2.1300952) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39078113) q[0];
sx q[0];
rz(-1.823714) q[0];
sx q[0];
rz(3.0506296) q[0];
rz(-2.2816351) q[1];
sx q[1];
rz(-2.0188315) q[1];
sx q[1];
rz(-1.8213173) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0866213) q[0];
sx q[0];
rz(-1.939659) q[0];
sx q[0];
rz(0.47199179) q[0];
x q[1];
rz(-0.57029057) q[2];
sx q[2];
rz(-1.4207134) q[2];
sx q[2];
rz(2.1652086) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7399866) q[1];
sx q[1];
rz(-1.1005797) q[1];
sx q[1];
rz(-0.92958881) q[1];
rz(-pi) q[2];
rz(0.28989132) q[3];
sx q[3];
rz(-0.93512669) q[3];
sx q[3];
rz(-1.0155201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0992574) q[2];
sx q[2];
rz(-2.1741185) q[2];
sx q[2];
rz(-0.60097224) q[2];
rz(2.6565334) q[3];
sx q[3];
rz(-0.22189134) q[3];
sx q[3];
rz(-1.4453567) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27959529) q[0];
sx q[0];
rz(-1.1789362) q[0];
sx q[0];
rz(-0.55554187) q[0];
rz(-0.034596054) q[1];
sx q[1];
rz(-2.3831773) q[1];
sx q[1];
rz(1.7506036) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9075508) q[0];
sx q[0];
rz(-2.4939135) q[0];
sx q[0];
rz(-0.56869047) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.820307) q[2];
sx q[2];
rz(-1.2148641) q[2];
sx q[2];
rz(-2.3552259) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78087378) q[1];
sx q[1];
rz(-2.2370173) q[1];
sx q[1];
rz(0.16155508) q[1];
x q[2];
rz(-1.6452541) q[3];
sx q[3];
rz(-0.7965318) q[3];
sx q[3];
rz(1.9394685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3283219) q[2];
sx q[2];
rz(-2.6997824) q[2];
sx q[2];
rz(0.39548809) q[2];
rz(1.8528806) q[3];
sx q[3];
rz(-1.6059395) q[3];
sx q[3];
rz(2.4718463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.678858) q[0];
sx q[0];
rz(-0.33927074) q[0];
sx q[0];
rz(-1.4920374) q[0];
rz(0.95343268) q[1];
sx q[1];
rz(-1.1089193) q[1];
sx q[1];
rz(-1.4377726) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2742845) q[0];
sx q[0];
rz(-2.0831997) q[0];
sx q[0];
rz(0.209765) q[0];
rz(0.17022325) q[2];
sx q[2];
rz(-1.7049689) q[2];
sx q[2];
rz(1.2302878) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7319665) q[1];
sx q[1];
rz(-2.1831174) q[1];
sx q[1];
rz(-0.97169505) q[1];
x q[2];
rz(3.1215454) q[3];
sx q[3];
rz(-1.123395) q[3];
sx q[3];
rz(-2.945154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64547223) q[2];
sx q[2];
rz(-2.7068553) q[2];
sx q[2];
rz(-2.9628741) q[2];
rz(-0.86137613) q[3];
sx q[3];
rz(-1.2025611) q[3];
sx q[3];
rz(-0.3716968) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9361967) q[0];
sx q[0];
rz(-1.9367171) q[0];
sx q[0];
rz(1.0937011) q[0];
rz(-2.4049092) q[1];
sx q[1];
rz(-1.271558) q[1];
sx q[1];
rz(-1.0587943) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98204389) q[0];
sx q[0];
rz(-1.0315572) q[0];
sx q[0];
rz(-0.3190785) q[0];
rz(-pi) q[1];
rz(2.2718614) q[2];
sx q[2];
rz(-0.5898925) q[2];
sx q[2];
rz(1.9260977) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0000856) q[1];
sx q[1];
rz(-2.7440789) q[1];
sx q[1];
rz(-1.8915528) q[1];
x q[2];
rz(-0.69370725) q[3];
sx q[3];
rz(-0.5872763) q[3];
sx q[3];
rz(1.8380084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.27292353) q[2];
sx q[2];
rz(-2.4464641) q[2];
sx q[2];
rz(1.6607364) q[2];
rz(2.7311834) q[3];
sx q[3];
rz(-1.4199665) q[3];
sx q[3];
rz(-2.9836392) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2720298) q[0];
sx q[0];
rz(-2.8700888) q[0];
sx q[0];
rz(-2.8503382) q[0];
rz(-0.60925305) q[1];
sx q[1];
rz(-1.6758502) q[1];
sx q[1];
rz(-1.7094918) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2925443) q[0];
sx q[0];
rz(-1.8869072) q[0];
sx q[0];
rz(-0.51878099) q[0];
x q[1];
rz(-2.6860793) q[2];
sx q[2];
rz(-1.3580772) q[2];
sx q[2];
rz(2.6754975) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6546302) q[1];
sx q[1];
rz(-1.527206) q[1];
sx q[1];
rz(-2.589588) q[1];
rz(-pi) q[2];
rz(-0.087177353) q[3];
sx q[3];
rz(-2.0012337) q[3];
sx q[3];
rz(-0.038539683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.46135205) q[2];
sx q[2];
rz(-2.0246918) q[2];
sx q[2];
rz(-2.0001901) q[2];
rz(-1.5348148) q[3];
sx q[3];
rz(-1.1780058) q[3];
sx q[3];
rz(0.46943584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5466945) q[0];
sx q[0];
rz(-1.6225157) q[0];
sx q[0];
rz(-1.7058104) q[0];
rz(-2.7728511) q[1];
sx q[1];
rz(-1.8992966) q[1];
sx q[1];
rz(-0.13175838) q[1];
rz(1.3651866) q[2];
sx q[2];
rz(-1.793135) q[2];
sx q[2];
rz(-1.8829913) q[2];
rz(1.7678123) q[3];
sx q[3];
rz(-1.1548629) q[3];
sx q[3];
rz(-2.9668273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
