OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7856287) q[0];
sx q[0];
rz(-0.2956737) q[0];
sx q[0];
rz(0.41391882) q[0];
rz(-2.8726481) q[1];
sx q[1];
rz(-2.7719331) q[1];
sx q[1];
rz(1.2764021) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1838754) q[0];
sx q[0];
rz(-1.4801637) q[0];
sx q[0];
rz(-3.0498226) q[0];
x q[1];
rz(-1.6501313) q[2];
sx q[2];
rz(-0.66986194) q[2];
sx q[2];
rz(-1.2487433) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2416083) q[1];
sx q[1];
rz(-1.2869724) q[1];
sx q[1];
rz(-2.349616) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5475575) q[3];
sx q[3];
rz(-1.7230936) q[3];
sx q[3];
rz(1.1588948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.74497491) q[2];
sx q[2];
rz(-1.2367542) q[2];
sx q[2];
rz(1.2006987) q[2];
rz(0.68136224) q[3];
sx q[3];
rz(-1.5228289) q[3];
sx q[3];
rz(-2.7549226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.50614) q[0];
sx q[0];
rz(-2.8400087) q[0];
sx q[0];
rz(-0.26148456) q[0];
rz(-1.6188949) q[1];
sx q[1];
rz(-2.6324184) q[1];
sx q[1];
rz(-1.8656628) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7000293) q[0];
sx q[0];
rz(-1.4744548) q[0];
sx q[0];
rz(-2.5896713) q[0];
rz(-pi) q[1];
rz(2.4894756) q[2];
sx q[2];
rz(-2.0154833) q[2];
sx q[2];
rz(2.7114078) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6258982) q[1];
sx q[1];
rz(-1.0401023) q[1];
sx q[1];
rz(1.9489392) q[1];
x q[2];
rz(-3.1072282) q[3];
sx q[3];
rz(-1.7196764) q[3];
sx q[3];
rz(-1.9158165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2168938) q[2];
sx q[2];
rz(-1.0328707) q[2];
sx q[2];
rz(0.72803289) q[2];
rz(1.3701471) q[3];
sx q[3];
rz(-2.5244505) q[3];
sx q[3];
rz(3.106015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.034123357) q[0];
sx q[0];
rz(-1.1119482) q[0];
sx q[0];
rz(0.6849826) q[0];
rz(-2.0383539) q[1];
sx q[1];
rz(-2.5220242) q[1];
sx q[1];
rz(1.0122976) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89268003) q[0];
sx q[0];
rz(-1.9143701) q[0];
sx q[0];
rz(1.1569886) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7660161) q[2];
sx q[2];
rz(-1.3998271) q[2];
sx q[2];
rz(-2.7107216) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5372181) q[1];
sx q[1];
rz(-0.24362016) q[1];
sx q[1];
rz(2.591406) q[1];
rz(-pi) q[2];
rz(-0.52148444) q[3];
sx q[3];
rz(-2.3064724) q[3];
sx q[3];
rz(-2.2948752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.95969069) q[2];
sx q[2];
rz(-1.5519374) q[2];
sx q[2];
rz(2.336179) q[2];
rz(2.3253333) q[3];
sx q[3];
rz(-2.5160242) q[3];
sx q[3];
rz(3.0474385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1212946) q[0];
sx q[0];
rz(-0.90885201) q[0];
sx q[0];
rz(0.62776172) q[0];
rz(-0.10920564) q[1];
sx q[1];
rz(-0.86995482) q[1];
sx q[1];
rz(-2.3039718) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8872466) q[0];
sx q[0];
rz(-1.5043048) q[0];
sx q[0];
rz(-2.6496135) q[0];
rz(-pi) q[1];
rz(-2.9170981) q[2];
sx q[2];
rz(-0.98735972) q[2];
sx q[2];
rz(-0.79085858) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.301103) q[1];
sx q[1];
rz(-0.87574848) q[1];
sx q[1];
rz(-1.3962505) q[1];
x q[2];
rz(-1.2628951) q[3];
sx q[3];
rz(-1.2114204) q[3];
sx q[3];
rz(1.9482159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8236905) q[2];
sx q[2];
rz(-2.6573942) q[2];
sx q[2];
rz(2.79706) q[2];
rz(-1.4065929) q[3];
sx q[3];
rz(-2.78648) q[3];
sx q[3];
rz(3.0936892) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.542881) q[0];
sx q[0];
rz(-2.5051835) q[0];
sx q[0];
rz(-2.6357546) q[0];
rz(1.0525616) q[1];
sx q[1];
rz(-2.274175) q[1];
sx q[1];
rz(-0.94598407) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1051798) q[0];
sx q[0];
rz(-2.7650021) q[0];
sx q[0];
rz(0.74684493) q[0];
x q[1];
rz(2.7840678) q[2];
sx q[2];
rz(-1.7058183) q[2];
sx q[2];
rz(0.46070489) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9717488) q[1];
sx q[1];
rz(-1.1976114) q[1];
sx q[1];
rz(-2.5118474) q[1];
rz(-pi) q[2];
rz(-3.0732156) q[3];
sx q[3];
rz(-1.42236) q[3];
sx q[3];
rz(2.1297586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0537009) q[2];
sx q[2];
rz(-2.044007) q[2];
sx q[2];
rz(1.5637448) q[2];
rz(1.0846694) q[3];
sx q[3];
rz(-2.2659149) q[3];
sx q[3];
rz(-2.5904371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56187335) q[0];
sx q[0];
rz(-2.5358574) q[0];
sx q[0];
rz(-2.9015923) q[0];
rz(0.9217681) q[1];
sx q[1];
rz(-2.0926937) q[1];
sx q[1];
rz(-1.9711432) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.474139) q[0];
sx q[0];
rz(-0.57155656) q[0];
sx q[0];
rz(3.0616687) q[0];
rz(-pi) q[1];
rz(0.35640772) q[2];
sx q[2];
rz(-1.4517759) q[2];
sx q[2];
rz(1.0072034) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.45615087) q[1];
sx q[1];
rz(-0.69760453) q[1];
sx q[1];
rz(2.3122723) q[1];
x q[2];
rz(1.7975054) q[3];
sx q[3];
rz(-0.41019687) q[3];
sx q[3];
rz(2.6997379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3971098) q[2];
sx q[2];
rz(-0.8195256) q[2];
sx q[2];
rz(2.5605555) q[2];
rz(2.2033612) q[3];
sx q[3];
rz(-0.56648985) q[3];
sx q[3];
rz(2.4375622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0051603) q[0];
sx q[0];
rz(-0.67661023) q[0];
sx q[0];
rz(-0.50049472) q[0];
rz(-2.9035134) q[1];
sx q[1];
rz(-1.4512738) q[1];
sx q[1];
rz(2.4949825) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2867271) q[0];
sx q[0];
rz(-1.5359794) q[0];
sx q[0];
rz(1.5725732) q[0];
rz(1.4562723) q[2];
sx q[2];
rz(-1.6045286) q[2];
sx q[2];
rz(-2.9303868) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2518066) q[1];
sx q[1];
rz(-1.4507035) q[1];
sx q[1];
rz(2.0577862) q[1];
rz(1.4073237) q[3];
sx q[3];
rz(-0.36296926) q[3];
sx q[3];
rz(0.27856871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.87767345) q[2];
sx q[2];
rz(-2.4223902) q[2];
sx q[2];
rz(2.4379697) q[2];
rz(-1.8536812) q[3];
sx q[3];
rz(-1.0602919) q[3];
sx q[3];
rz(-2.6357486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8877761) q[0];
sx q[0];
rz(-1.1340589) q[0];
sx q[0];
rz(1.3251086) q[0];
rz(-0.77249402) q[1];
sx q[1];
rz(-2.019181) q[1];
sx q[1];
rz(0.17556369) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52161509) q[0];
sx q[0];
rz(-2.3562263) q[0];
sx q[0];
rz(0.26544858) q[0];
rz(-pi) q[1];
rz(-0.63736659) q[2];
sx q[2];
rz(-0.91203472) q[2];
sx q[2];
rz(0.8085685) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0766636) q[1];
sx q[1];
rz(-2.6573219) q[1];
sx q[1];
rz(1.8736429) q[1];
rz(2.9977071) q[3];
sx q[3];
rz(-0.7215313) q[3];
sx q[3];
rz(-0.084621457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0385711) q[2];
sx q[2];
rz(-1.0744289) q[2];
sx q[2];
rz(-2.0477022) q[2];
rz(1.966018) q[3];
sx q[3];
rz(-0.75785494) q[3];
sx q[3];
rz(0.28483835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4232101) q[0];
sx q[0];
rz(-0.045561401) q[0];
sx q[0];
rz(-1.3257931) q[0];
rz(2.6153053) q[1];
sx q[1];
rz(-2.278639) q[1];
sx q[1];
rz(-2.3525499) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1068971) q[0];
sx q[0];
rz(-1.7321087) q[0];
sx q[0];
rz(3.0437896) q[0];
rz(-pi) q[1];
rz(1.3157719) q[2];
sx q[2];
rz(-1.2256283) q[2];
sx q[2];
rz(1.5288625) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.55780667) q[1];
sx q[1];
rz(-2.4135488) q[1];
sx q[1];
rz(-2.1017722) q[1];
rz(-2.2317722) q[3];
sx q[3];
rz(-2.2299521) q[3];
sx q[3];
rz(-0.59162635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9290756) q[2];
sx q[2];
rz(-2.2079461) q[2];
sx q[2];
rz(0.90544256) q[2];
rz(-2.2255911) q[3];
sx q[3];
rz(-1.6986877) q[3];
sx q[3];
rz(-1.0441095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4395897) q[0];
sx q[0];
rz(-2.3083394) q[0];
sx q[0];
rz(-0.92779094) q[0];
rz(1.0935498) q[1];
sx q[1];
rz(-2.5826726) q[1];
sx q[1];
rz(3.0082026) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095208406) q[0];
sx q[0];
rz(-2.781087) q[0];
sx q[0];
rz(0.62490873) q[0];
rz(0.9993365) q[2];
sx q[2];
rz(-1.7448057) q[2];
sx q[2];
rz(-2.0561754) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.16001734) q[1];
sx q[1];
rz(-2.3381553) q[1];
sx q[1];
rz(2.395242) q[1];
rz(2.0285151) q[3];
sx q[3];
rz(-1.3821756) q[3];
sx q[3];
rz(2.6831804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5043958) q[2];
sx q[2];
rz(-1.7375676) q[2];
sx q[2];
rz(1.8782328) q[2];
rz(1.0697621) q[3];
sx q[3];
rz(-2.1061149) q[3];
sx q[3];
rz(3.0647762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1509811) q[0];
sx q[0];
rz(-2.4438416) q[0];
sx q[0];
rz(-1.7846815) q[0];
rz(-1.4487343) q[1];
sx q[1];
rz(-0.81304638) q[1];
sx q[1];
rz(-0.1437694) q[1];
rz(0.31126304) q[2];
sx q[2];
rz(-1.0672206) q[2];
sx q[2];
rz(-0.5819566) q[2];
rz(0.24604194) q[3];
sx q[3];
rz(-1.9046219) q[3];
sx q[3];
rz(0.73840284) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
