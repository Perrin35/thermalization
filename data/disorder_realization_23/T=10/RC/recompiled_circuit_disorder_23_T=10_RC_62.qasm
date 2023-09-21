OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2620579) q[0];
sx q[0];
rz(-1.7320002) q[0];
sx q[0];
rz(1.4341266) q[0];
rz(0.6342451) q[1];
sx q[1];
rz(-2.5399962) q[1];
sx q[1];
rz(-0.4184202) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.943676) q[0];
sx q[0];
rz(-1.1063873) q[0];
sx q[0];
rz(0.29505131) q[0];
rz(-pi) q[1];
rz(-3.0298972) q[2];
sx q[2];
rz(-1.7061491) q[2];
sx q[2];
rz(0.40458194) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.99416713) q[1];
sx q[1];
rz(-0.52417437) q[1];
sx q[1];
rz(1.0679354) q[1];
x q[2];
rz(2.9207346) q[3];
sx q[3];
rz(-2.0005895) q[3];
sx q[3];
rz(2.5151099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3216386) q[2];
sx q[2];
rz(-1.3691838) q[2];
sx q[2];
rz(2.3036172) q[2];
rz(-0.49301246) q[3];
sx q[3];
rz(-2.8686782) q[3];
sx q[3];
rz(-3.0626007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3943966) q[0];
sx q[0];
rz(-2.4173739) q[0];
sx q[0];
rz(1.2778506) q[0];
rz(0.17678075) q[1];
sx q[1];
rz(-1.3143833) q[1];
sx q[1];
rz(-0.4321672) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9958187) q[0];
sx q[0];
rz(-2.5403025) q[0];
sx q[0];
rz(0.50253089) q[0];
rz(-1.1439267) q[2];
sx q[2];
rz(-1.9383213) q[2];
sx q[2];
rz(-1.9232242) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.13465263) q[1];
sx q[1];
rz(-1.4538527) q[1];
sx q[1];
rz(-2.0140531) q[1];
x q[2];
rz(0.024859419) q[3];
sx q[3];
rz(-1.0059352) q[3];
sx q[3];
rz(-2.3853175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.54923487) q[2];
sx q[2];
rz(-1.2640415) q[2];
sx q[2];
rz(-2.6548927) q[2];
rz(-1.7633847) q[3];
sx q[3];
rz(-1.2599726) q[3];
sx q[3];
rz(2.6087705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.10087) q[0];
sx q[0];
rz(-2.3951055) q[0];
sx q[0];
rz(-2.7242463) q[0];
rz(-1.6529282) q[1];
sx q[1];
rz(-0.54549837) q[1];
sx q[1];
rz(2.6352077) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9837512) q[0];
sx q[0];
rz(-2.7451773) q[0];
sx q[0];
rz(-2.0435964) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8110397) q[2];
sx q[2];
rz(-2.2908387) q[2];
sx q[2];
rz(-0.70149295) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3370034) q[1];
sx q[1];
rz(-0.5491921) q[1];
sx q[1];
rz(-2.5227491) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9403463) q[3];
sx q[3];
rz(-1.7412211) q[3];
sx q[3];
rz(2.2263262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.69904077) q[2];
sx q[2];
rz(-0.46135819) q[2];
sx q[2];
rz(-0.59147269) q[2];
rz(0.58602035) q[3];
sx q[3];
rz(-1.932671) q[3];
sx q[3];
rz(1.4311786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1699003) q[0];
sx q[0];
rz(-2.5132892) q[0];
sx q[0];
rz(-0.83918321) q[0];
rz(-3.1160141) q[1];
sx q[1];
rz(-0.69568101) q[1];
sx q[1];
rz(1.5485839) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4913113) q[0];
sx q[0];
rz(-2.6822753) q[0];
sx q[0];
rz(-1.7772872) q[0];
rz(0.71009212) q[2];
sx q[2];
rz(-0.93821628) q[2];
sx q[2];
rz(-1.8284423) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3241987) q[1];
sx q[1];
rz(-0.84016582) q[1];
sx q[1];
rz(-0.52449951) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5769916) q[3];
sx q[3];
rz(-0.49800107) q[3];
sx q[3];
rz(-2.6548487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8362391) q[2];
sx q[2];
rz(-2.2576015) q[2];
sx q[2];
rz(0.099686064) q[2];
rz(-0.95885197) q[3];
sx q[3];
rz(-1.3189664) q[3];
sx q[3];
rz(-1.8166186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5754159) q[0];
sx q[0];
rz(-1.382099) q[0];
sx q[0];
rz(0.25594041) q[0];
rz(-0.4610962) q[1];
sx q[1];
rz(-1.0436811) q[1];
sx q[1];
rz(-2.3815313) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5417267) q[0];
sx q[0];
rz(-1.6875629) q[0];
sx q[0];
rz(1.4670502) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9854457) q[2];
sx q[2];
rz(-1.5125325) q[2];
sx q[2];
rz(-2.0410048) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.252994) q[1];
sx q[1];
rz(-2.9610486) q[1];
sx q[1];
rz(-2.351159) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0403967) q[3];
sx q[3];
rz(-0.70221516) q[3];
sx q[3];
rz(-1.5576253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.57050675) q[2];
sx q[2];
rz(-1.7692302) q[2];
sx q[2];
rz(-2.467353) q[2];
rz(2.9267866) q[3];
sx q[3];
rz(-0.45682296) q[3];
sx q[3];
rz(0.017344346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.53133416) q[0];
sx q[0];
rz(-1.6711618) q[0];
sx q[0];
rz(-1.9624788) q[0];
rz(-0.20482652) q[1];
sx q[1];
rz(-2.3463459) q[1];
sx q[1];
rz(-2.0746322) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93766312) q[0];
sx q[0];
rz(-1.5563037) q[0];
sx q[0];
rz(0.020676215) q[0];
rz(-pi) q[1];
rz(-0.80337556) q[2];
sx q[2];
rz(-0.57758812) q[2];
sx q[2];
rz(-0.20882777) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8384335) q[1];
sx q[1];
rz(-1.618297) q[1];
sx q[1];
rz(2.3179503) q[1];
rz(1.7339891) q[3];
sx q[3];
rz(-0.7810775) q[3];
sx q[3];
rz(-1.2490602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.71010464) q[2];
sx q[2];
rz(-1.8523214) q[2];
sx q[2];
rz(-2.5816494) q[2];
rz(-2.4152749) q[3];
sx q[3];
rz(-2.8328219) q[3];
sx q[3];
rz(2.8360951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85754919) q[0];
sx q[0];
rz(-2.5950268) q[0];
sx q[0];
rz(1.7204826) q[0];
rz(-0.20206085) q[1];
sx q[1];
rz(-1.4338564) q[1];
sx q[1];
rz(0.85817671) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1976397) q[0];
sx q[0];
rz(-1.7174935) q[0];
sx q[0];
rz(3.0118045) q[0];
rz(-pi) q[1];
rz(-1.6740587) q[2];
sx q[2];
rz(-0.41245663) q[2];
sx q[2];
rz(2.7445284) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8691751) q[1];
sx q[1];
rz(-0.032748001) q[1];
sx q[1];
rz(-0.48519022) q[1];
x q[2];
rz(0.51123294) q[3];
sx q[3];
rz(-0.63496642) q[3];
sx q[3];
rz(-2.883203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.28785607) q[2];
sx q[2];
rz(-2.6617472) q[2];
sx q[2];
rz(-1.8161592) q[2];
rz(2.251513) q[3];
sx q[3];
rz(-1.9944913) q[3];
sx q[3];
rz(0.98852283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4637852) q[0];
sx q[0];
rz(-0.57254922) q[0];
sx q[0];
rz(-2.7668787) q[0];
rz(0.97887865) q[1];
sx q[1];
rz(-2.4596877) q[1];
sx q[1];
rz(1.3495061) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23096202) q[0];
sx q[0];
rz(-0.9696784) q[0];
sx q[0];
rz(-2.0448951) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1153203) q[2];
sx q[2];
rz(-0.71825829) q[2];
sx q[2];
rz(2.9422613) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0950349) q[1];
sx q[1];
rz(-1.7525502) q[1];
sx q[1];
rz(-0.05038105) q[1];
rz(-pi) q[2];
rz(-2.5958063) q[3];
sx q[3];
rz(-1.3457527) q[3];
sx q[3];
rz(-0.033566098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.65138856) q[2];
sx q[2];
rz(-2.3888402) q[2];
sx q[2];
rz(-2.6728969) q[2];
rz(-1.1941341) q[3];
sx q[3];
rz(-1.3041376) q[3];
sx q[3];
rz(-1.3635925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63012183) q[0];
sx q[0];
rz(-0.90634316) q[0];
sx q[0];
rz(-1.2619031) q[0];
rz(0.17503861) q[1];
sx q[1];
rz(-1.9997528) q[1];
sx q[1];
rz(1.6040241) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3106829) q[0];
sx q[0];
rz(-1.791782) q[0];
sx q[0];
rz(-1.2587147) q[0];
rz(-pi) q[1];
rz(2.9109355) q[2];
sx q[2];
rz(-0.42867491) q[2];
sx q[2];
rz(-2.8187657) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9201587) q[1];
sx q[1];
rz(-2.4281574) q[1];
sx q[1];
rz(-0.17381298) q[1];
rz(-pi) q[2];
rz(-0.92026199) q[3];
sx q[3];
rz(-2.5839621) q[3];
sx q[3];
rz(-2.695431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1016772) q[2];
sx q[2];
rz(-0.95255178) q[2];
sx q[2];
rz(0.76134479) q[2];
rz(-2.2411761) q[3];
sx q[3];
rz(-2.5420928) q[3];
sx q[3];
rz(3.0925687) q[3];
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
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.802357) q[0];
sx q[0];
rz(-0.46827066) q[0];
sx q[0];
rz(0.21690579) q[0];
rz(0.63198173) q[1];
sx q[1];
rz(-1.6501553) q[1];
sx q[1];
rz(2.1868618) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5548984) q[0];
sx q[0];
rz(-2.2838755) q[0];
sx q[0];
rz(-1.1416392) q[0];
x q[1];
rz(-1.9913313) q[2];
sx q[2];
rz(-1.0610126) q[2];
sx q[2];
rz(-1.860294) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54143822) q[1];
sx q[1];
rz(-2.6791875) q[1];
sx q[1];
rz(2.9701783) q[1];
x q[2];
rz(-2.2050489) q[3];
sx q[3];
rz(-1.9359971) q[3];
sx q[3];
rz(1.6983502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1251936) q[2];
sx q[2];
rz(-1.23896) q[2];
sx q[2];
rz(2.1968502) q[2];
rz(0.38481209) q[3];
sx q[3];
rz(-2.0231569) q[3];
sx q[3];
rz(2.184536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.64086296) q[0];
sx q[0];
rz(-2.6364115) q[0];
sx q[0];
rz(-1.5873948) q[0];
rz(-2.2568933) q[1];
sx q[1];
rz(-0.90507602) q[1];
sx q[1];
rz(-0.25837635) q[1];
rz(2.2434071) q[2];
sx q[2];
rz(-0.85815103) q[2];
sx q[2];
rz(1.382538) q[2];
rz(-0.20529071) q[3];
sx q[3];
rz(-1.1355573) q[3];
sx q[3];
rz(-0.60187403) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];