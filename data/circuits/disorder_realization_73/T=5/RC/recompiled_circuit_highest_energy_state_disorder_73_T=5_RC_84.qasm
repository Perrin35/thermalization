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
rz(-2.7849164) q[0];
sx q[0];
rz(-1.4465605) q[0];
sx q[0];
rz(0.90176982) q[0];
rz(4.4737368) q[1];
sx q[1];
rz(2.6982215) q[1];
sx q[1];
rz(8.659521) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1674102) q[0];
sx q[0];
rz(-1.2077792) q[0];
sx q[0];
rz(-0.98591759) q[0];
rz(-pi) q[1];
rz(1.9087431) q[2];
sx q[2];
rz(-1.3409919) q[2];
sx q[2];
rz(2.1462962) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3292685) q[1];
sx q[1];
rz(-1.6109038) q[1];
sx q[1];
rz(0.33848156) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0115143) q[3];
sx q[3];
rz(-1.2941232) q[3];
sx q[3];
rz(-2.193146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.29609933) q[2];
sx q[2];
rz(-2.0534434) q[2];
sx q[2];
rz(-1.4720526) q[2];
rz(-1.7146401) q[3];
sx q[3];
rz(-1.497437) q[3];
sx q[3];
rz(0.53372395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.958309) q[0];
sx q[0];
rz(-1.6338209) q[0];
sx q[0];
rz(-0.84877745) q[0];
rz(-0.035004184) q[1];
sx q[1];
rz(-1.1468381) q[1];
sx q[1];
rz(2.1880207) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35915056) q[0];
sx q[0];
rz(-2.0277362) q[0];
sx q[0];
rz(0.63553973) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9601279) q[2];
sx q[2];
rz(-1.2876533) q[2];
sx q[2];
rz(2.2623537) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2767859) q[1];
sx q[1];
rz(-1.3753483) q[1];
sx q[1];
rz(-1.2246183) q[1];
rz(-pi) q[2];
rz(2.4783432) q[3];
sx q[3];
rz(-0.8732855) q[3];
sx q[3];
rz(2.8467941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1635052) q[2];
sx q[2];
rz(-1.7814025) q[2];
sx q[2];
rz(-2.2920442) q[2];
rz(-2.9761525) q[3];
sx q[3];
rz(-1.6980419) q[3];
sx q[3];
rz(-2.3918242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7716832) q[0];
sx q[0];
rz(-2.219438) q[0];
sx q[0];
rz(-2.7396696) q[0];
rz(-2.9882714) q[1];
sx q[1];
rz(-0.17686495) q[1];
sx q[1];
rz(2.0053999) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76321852) q[0];
sx q[0];
rz(-2.741886) q[0];
sx q[0];
rz(1.6016774) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6499475) q[2];
sx q[2];
rz(-0.53722135) q[2];
sx q[2];
rz(1.1342837) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2991645) q[1];
sx q[1];
rz(-1.8106726) q[1];
sx q[1];
rz(0.23934083) q[1];
rz(-pi) q[2];
rz(2.4506684) q[3];
sx q[3];
rz(-2.5657585) q[3];
sx q[3];
rz(1.761328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.084126964) q[2];
sx q[2];
rz(-0.87856138) q[2];
sx q[2];
rz(-3.086536) q[2];
rz(1.7940686) q[3];
sx q[3];
rz(-1.3301347) q[3];
sx q[3];
rz(-2.1373035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36443001) q[0];
sx q[0];
rz(-0.20474064) q[0];
sx q[0];
rz(-1.5740016) q[0];
rz(-2.4500627) q[1];
sx q[1];
rz(-2.4722996) q[1];
sx q[1];
rz(2.9561668) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2381993) q[0];
sx q[0];
rz(-1.1739587) q[0];
sx q[0];
rz(1.7522041) q[0];
rz(-2.4302684) q[2];
sx q[2];
rz(-0.60718107) q[2];
sx q[2];
rz(-1.173623) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.54900995) q[1];
sx q[1];
rz(-1.8235778) q[1];
sx q[1];
rz(1.8922674) q[1];
x q[2];
rz(-2.809285) q[3];
sx q[3];
rz(-0.57113591) q[3];
sx q[3];
rz(0.50779282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4039679) q[2];
sx q[2];
rz(-1.4289958) q[2];
sx q[2];
rz(-3.1363623) q[2];
rz(2.7541006) q[3];
sx q[3];
rz(-0.5046851) q[3];
sx q[3];
rz(-1.8385878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72652793) q[0];
sx q[0];
rz(-0.98353493) q[0];
sx q[0];
rz(2.5256185) q[0];
rz(-1.9527324) q[1];
sx q[1];
rz(-2.523246) q[1];
sx q[1];
rz(-2.628285) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6270638) q[0];
sx q[0];
rz(-1.5634057) q[0];
sx q[0];
rz(-1.2891931) q[0];
x q[1];
rz(-1.2890069) q[2];
sx q[2];
rz(-0.86494397) q[2];
sx q[2];
rz(-0.47823804) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3860046) q[1];
sx q[1];
rz(-1.6428448) q[1];
sx q[1];
rz(2.5064038) q[1];
rz(2.0381472) q[3];
sx q[3];
rz(-1.6765521) q[3];
sx q[3];
rz(-2.8018564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1335699) q[2];
sx q[2];
rz(-1.8654537) q[2];
sx q[2];
rz(-2.3991154) q[2];
rz(1.4437458) q[3];
sx q[3];
rz(-1.7323114) q[3];
sx q[3];
rz(1.1013364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3105069) q[0];
sx q[0];
rz(-2.0948912) q[0];
sx q[0];
rz(-2.781784) q[0];
rz(-2.2911435) q[1];
sx q[1];
rz(-1.9603739) q[1];
sx q[1];
rz(-2.6757619) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3208) q[0];
sx q[0];
rz(-0.75220901) q[0];
sx q[0];
rz(-1.2515995) q[0];
x q[1];
rz(-2.3783422) q[2];
sx q[2];
rz(-1.4723198) q[2];
sx q[2];
rz(1.5710448) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9656495) q[1];
sx q[1];
rz(-1.6217983) q[1];
sx q[1];
rz(-0.24105625) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96420607) q[3];
sx q[3];
rz(-1.0830942) q[3];
sx q[3];
rz(1.694198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3992074) q[2];
sx q[2];
rz(-1.9611605) q[2];
sx q[2];
rz(-0.36435374) q[2];
rz(-2.897701) q[3];
sx q[3];
rz(-3.0920005) q[3];
sx q[3];
rz(1.9243139) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57109443) q[0];
sx q[0];
rz(-0.97935337) q[0];
sx q[0];
rz(-1.97557) q[0];
rz(2.1265538) q[1];
sx q[1];
rz(-2.6045585) q[1];
sx q[1];
rz(-2.7551415) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58182654) q[0];
sx q[0];
rz(-0.21187267) q[0];
sx q[0];
rz(-1.7170402) q[0];
x q[1];
rz(-2.831424) q[2];
sx q[2];
rz(-1.3259238) q[2];
sx q[2];
rz(0.32178165) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9484752) q[1];
sx q[1];
rz(-1.5962692) q[1];
sx q[1];
rz(1.6465882) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0160955) q[3];
sx q[3];
rz(-1.0319796) q[3];
sx q[3];
rz(-1.5244689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.40547907) q[2];
sx q[2];
rz(-1.4961286) q[2];
sx q[2];
rz(0.22769895) q[2];
rz(2.0685711) q[3];
sx q[3];
rz(-2.3927972) q[3];
sx q[3];
rz(-0.29357114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88184083) q[0];
sx q[0];
rz(-0.89983928) q[0];
sx q[0];
rz(-0.42801273) q[0];
rz(-2.1233066) q[1];
sx q[1];
rz(-2.7052453) q[1];
sx q[1];
rz(0.87103081) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6414531) q[0];
sx q[0];
rz(-0.79384168) q[0];
sx q[0];
rz(2.3026233) q[0];
x q[1];
rz(-2.9350014) q[2];
sx q[2];
rz(-1.2204514) q[2];
sx q[2];
rz(-1.3526582) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.97058875) q[1];
sx q[1];
rz(-1.4297545) q[1];
sx q[1];
rz(-3.019437) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4860367) q[3];
sx q[3];
rz(-2.3456315) q[3];
sx q[3];
rz(-1.9565005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.009306) q[2];
sx q[2];
rz(-0.79557482) q[2];
sx q[2];
rz(-1.2786678) q[2];
rz(2.5631185) q[3];
sx q[3];
rz(-0.63251248) q[3];
sx q[3];
rz(1.1540029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4671675) q[0];
sx q[0];
rz(-2.9245057) q[0];
sx q[0];
rz(0.59980741) q[0];
rz(-1.2990052) q[1];
sx q[1];
rz(-1.5312342) q[1];
sx q[1];
rz(-0.64186796) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3369747) q[0];
sx q[0];
rz(-2.0463159) q[0];
sx q[0];
rz(1.9561951) q[0];
rz(2.9974143) q[2];
sx q[2];
rz(-0.81071172) q[2];
sx q[2];
rz(1.7264896) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7781592) q[1];
sx q[1];
rz(-1.5892519) q[1];
sx q[1];
rz(-0.70779558) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3849726) q[3];
sx q[3];
rz(-1.1706268) q[3];
sx q[3];
rz(0.64368806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.32144) q[2];
sx q[2];
rz(-1.4250616) q[2];
sx q[2];
rz(-2.5313306) q[2];
rz(2.6214456) q[3];
sx q[3];
rz(-2.0870233) q[3];
sx q[3];
rz(2.648073) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.801012) q[0];
sx q[0];
rz(-1.8980674) q[0];
sx q[0];
rz(-0.93801671) q[0];
rz(2.7998789) q[1];
sx q[1];
rz(-1.7594756) q[1];
sx q[1];
rz(2.9843073) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.324824) q[0];
sx q[0];
rz(-0.95476645) q[0];
sx q[0];
rz(-1.8094713) q[0];
rz(-0.75802676) q[2];
sx q[2];
rz(-1.5093414) q[2];
sx q[2];
rz(3.0830992) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.966519) q[1];
sx q[1];
rz(-2.0957114) q[1];
sx q[1];
rz(-1.2205842) q[1];
rz(-pi) q[2];
rz(0.25453849) q[3];
sx q[3];
rz(-1.1738281) q[3];
sx q[3];
rz(-2.1074668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9226795) q[2];
sx q[2];
rz(-0.34505406) q[2];
sx q[2];
rz(-1.867713) q[2];
rz(0.0056565469) q[3];
sx q[3];
rz(-1.4414682) q[3];
sx q[3];
rz(-2.1605632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8511178) q[0];
sx q[0];
rz(-1.1707476) q[0];
sx q[0];
rz(0.049402417) q[0];
rz(-2.6750917) q[1];
sx q[1];
rz(-1.3460881) q[1];
sx q[1];
rz(0.64429611) q[1];
rz(2.2412655) q[2];
sx q[2];
rz(-1.1569958) q[2];
sx q[2];
rz(-0.3668084) q[2];
rz(1.4884819) q[3];
sx q[3];
rz(-2.001279) q[3];
sx q[3];
rz(1.6829987) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
