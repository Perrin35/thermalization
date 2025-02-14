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
rz(1.264313) q[0];
sx q[0];
rz(6.6067196) q[0];
sx q[0];
rz(8.9759448) q[0];
rz(1.9857061) q[1];
sx q[1];
rz(-1.356025) q[1];
sx q[1];
rz(1.1745656) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33307459) q[0];
sx q[0];
rz(-1.9155598) q[0];
sx q[0];
rz(0.60057849) q[0];
x q[1];
rz(-1.5023793) q[2];
sx q[2];
rz(-1.2478154) q[2];
sx q[2];
rz(-3.0634865) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.9372896) q[1];
sx q[1];
rz(-2.1912592) q[1];
sx q[1];
rz(-0.043067769) q[1];
x q[2];
rz(0.17711344) q[3];
sx q[3];
rz(-2.6209954) q[3];
sx q[3];
rz(2.2442852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.590098) q[2];
sx q[2];
rz(-1.3000725) q[2];
sx q[2];
rz(0.91280118) q[2];
rz(1.270795) q[3];
sx q[3];
rz(-1.0400892) q[3];
sx q[3];
rz(-2.2916268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0047282334) q[0];
sx q[0];
rz(-2.5623463) q[0];
sx q[0];
rz(-2.7194523) q[0];
rz(0.36960754) q[1];
sx q[1];
rz(-2.044675) q[1];
sx q[1];
rz(-2.2014528) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63623171) q[0];
sx q[0];
rz(-0.32948895) q[0];
sx q[0];
rz(-1.5630653) q[0];
x q[1];
rz(1.3708417) q[2];
sx q[2];
rz(-2.5995289) q[2];
sx q[2];
rz(2.0393537) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89366363) q[1];
sx q[1];
rz(-1.0808696) q[1];
sx q[1];
rz(-0.72743261) q[1];
rz(-0.27249725) q[3];
sx q[3];
rz(-2.2292308) q[3];
sx q[3];
rz(2.8373425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9422354) q[2];
sx q[2];
rz(-2.4713559) q[2];
sx q[2];
rz(-0.90636903) q[2];
rz(0.97936112) q[3];
sx q[3];
rz(-2.7370079) q[3];
sx q[3];
rz(1.8766859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1449428) q[0];
sx q[0];
rz(-0.066983797) q[0];
sx q[0];
rz(-1.8970733) q[0];
rz(2.7527346) q[1];
sx q[1];
rz(-1.2797979) q[1];
sx q[1];
rz(-2.8188474) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9661315) q[0];
sx q[0];
rz(-1.4585988) q[0];
sx q[0];
rz(0.66288046) q[0];
x q[1];
rz(0.30827807) q[2];
sx q[2];
rz(-2.1090504) q[2];
sx q[2];
rz(-0.27632144) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6420329) q[1];
sx q[1];
rz(-0.82553464) q[1];
sx q[1];
rz(0.31742974) q[1];
x q[2];
rz(1.6204319) q[3];
sx q[3];
rz(-1.2200886) q[3];
sx q[3];
rz(-1.8010822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2286223) q[2];
sx q[2];
rz(-1.3801489) q[2];
sx q[2];
rz(-1.3529533) q[2];
rz(2.8690423) q[3];
sx q[3];
rz(-1.8504146) q[3];
sx q[3];
rz(1.6779617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6005741) q[0];
sx q[0];
rz(-2.8193642) q[0];
sx q[0];
rz(-1.898265) q[0];
rz(0.86878949) q[1];
sx q[1];
rz(-2.181535) q[1];
sx q[1];
rz(2.0713846) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3995099) q[0];
sx q[0];
rz(-0.94219452) q[0];
sx q[0];
rz(2.5058772) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44311503) q[2];
sx q[2];
rz(-0.84544467) q[2];
sx q[2];
rz(-2.1381174) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94302109) q[1];
sx q[1];
rz(-2.7987438) q[1];
sx q[1];
rz(-1.917761) q[1];
x q[2];
rz(1.7620874) q[3];
sx q[3];
rz(-1.9076348) q[3];
sx q[3];
rz(1.2113406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5493245) q[2];
sx q[2];
rz(-2.3531239) q[2];
sx q[2];
rz(0.92174021) q[2];
rz(-2.8978469) q[3];
sx q[3];
rz(-1.4344401) q[3];
sx q[3];
rz(2.8488979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.155596) q[0];
sx q[0];
rz(-1.8159741) q[0];
sx q[0];
rz(0.50088125) q[0];
rz(-1.1174551) q[1];
sx q[1];
rz(-1.3331579) q[1];
sx q[1];
rz(0.31455988) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9524744) q[0];
sx q[0];
rz(-0.50273147) q[0];
sx q[0];
rz(-2.377654) q[0];
rz(-pi) q[1];
rz(-1.4681513) q[2];
sx q[2];
rz(-0.043641239) q[2];
sx q[2];
rz(-2.1473644) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5774975) q[1];
sx q[1];
rz(-2.4397813) q[1];
sx q[1];
rz(-2.2086269) q[1];
x q[2];
rz(0.1184095) q[3];
sx q[3];
rz(-1.1573912) q[3];
sx q[3];
rz(-1.4865231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.306281) q[2];
sx q[2];
rz(-1.1113144) q[2];
sx q[2];
rz(-1.8446946) q[2];
rz(0.48526192) q[3];
sx q[3];
rz(-1.5596215) q[3];
sx q[3];
rz(2.0280793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6473963) q[0];
sx q[0];
rz(-1.8504471) q[0];
sx q[0];
rz(-0.099763481) q[0];
rz(-1.8708723) q[1];
sx q[1];
rz(-2.2427509) q[1];
sx q[1];
rz(-1.7521923) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6806599) q[0];
sx q[0];
rz(-1.4062728) q[0];
sx q[0];
rz(-0.27030944) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4657137) q[2];
sx q[2];
rz(-1.7438403) q[2];
sx q[2];
rz(1.4937173) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2459316) q[1];
sx q[1];
rz(-2.1626923) q[1];
sx q[1];
rz(0.81879692) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0148195) q[3];
sx q[3];
rz(-2.354458) q[3];
sx q[3];
rz(-1.8925557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.49326593) q[2];
sx q[2];
rz(-0.47097012) q[2];
sx q[2];
rz(-0.89034447) q[2];
rz(1.9503615) q[3];
sx q[3];
rz(-1.8343364) q[3];
sx q[3];
rz(0.90203917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7118199) q[0];
sx q[0];
rz(-3.0968102) q[0];
sx q[0];
rz(3.1088767) q[0];
rz(0.50321594) q[1];
sx q[1];
rz(-2.0423753) q[1];
sx q[1];
rz(3.018766) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3204725) q[0];
sx q[0];
rz(-0.80424612) q[0];
sx q[0];
rz(-0.45566092) q[0];
rz(-pi) q[1];
rz(-0.57374222) q[2];
sx q[2];
rz(-2.5649568) q[2];
sx q[2];
rz(-1.5337616) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1122857) q[1];
sx q[1];
rz(-2.0404792) q[1];
sx q[1];
rz(-0.11135688) q[1];
rz(-pi) q[2];
rz(2.7890754) q[3];
sx q[3];
rz(-1.4488701) q[3];
sx q[3];
rz(-0.33559775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.34933019) q[2];
sx q[2];
rz(-0.96668875) q[2];
sx q[2];
rz(-1.1518504) q[2];
rz(-0.88519111) q[3];
sx q[3];
rz(-1.3693634) q[3];
sx q[3];
rz(1.4256029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18871466) q[0];
sx q[0];
rz(-0.9335683) q[0];
sx q[0];
rz(-0.92920148) q[0];
rz(-0.62905351) q[1];
sx q[1];
rz(-1.355143) q[1];
sx q[1];
rz(-2.3984875) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4152195) q[0];
sx q[0];
rz(-2.6872244) q[0];
sx q[0];
rz(-1.3903862) q[0];
rz(-0.17878187) q[2];
sx q[2];
rz(-2.6465073) q[2];
sx q[2];
rz(2.924233) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8344017) q[1];
sx q[1];
rz(-1.9842923) q[1];
sx q[1];
rz(-1.6430699) q[1];
rz(-2.6499641) q[3];
sx q[3];
rz(-2.1475882) q[3];
sx q[3];
rz(2.4099809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5653845) q[2];
sx q[2];
rz(-2.3614063) q[2];
sx q[2];
rz(0.72594491) q[2];
rz(-2.7706326) q[3];
sx q[3];
rz(-1.5981263) q[3];
sx q[3];
rz(-2.7788739) q[3];
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
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87042701) q[0];
sx q[0];
rz(-2.3112264) q[0];
sx q[0];
rz(2.5776432) q[0];
rz(1.14934) q[1];
sx q[1];
rz(-2.1795858) q[1];
sx q[1];
rz(-0.74251485) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19631736) q[0];
sx q[0];
rz(-0.68500297) q[0];
sx q[0];
rz(2.3514868) q[0];
rz(-pi) q[1];
rz(-2.5385802) q[2];
sx q[2];
rz(-2.6508287) q[2];
sx q[2];
rz(3.1163851) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92723211) q[1];
sx q[1];
rz(-2.7013268) q[1];
sx q[1];
rz(3.0378067) q[1];
rz(-2.2015195) q[3];
sx q[3];
rz(-0.82118644) q[3];
sx q[3];
rz(1.8668777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5138862) q[2];
sx q[2];
rz(-1.0739001) q[2];
sx q[2];
rz(2.7970496) q[2];
rz(-2.6978317) q[3];
sx q[3];
rz(-1.0756451) q[3];
sx q[3];
rz(-1.8189323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3798856) q[0];
sx q[0];
rz(-0.2825309) q[0];
sx q[0];
rz(-0.66201061) q[0];
rz(-0.33540353) q[1];
sx q[1];
rz(-1.6796203) q[1];
sx q[1];
rz(-0.9368771) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59360817) q[0];
sx q[0];
rz(-2.4146705) q[0];
sx q[0];
rz(2.1288104) q[0];
x q[1];
rz(3.1179948) q[2];
sx q[2];
rz(-1.840971) q[2];
sx q[2];
rz(-0.23575704) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.42439991) q[1];
sx q[1];
rz(-0.48312995) q[1];
sx q[1];
rz(-2.4880954) q[1];
rz(-2.2912737) q[3];
sx q[3];
rz(-2.8788044) q[3];
sx q[3];
rz(-2.9966054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.718049) q[2];
sx q[2];
rz(-0.23423883) q[2];
sx q[2];
rz(-0.49918175) q[2];
rz(-1.4623803) q[3];
sx q[3];
rz(-1.1464109) q[3];
sx q[3];
rz(-1.4596938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41596169) q[0];
sx q[0];
rz(-1.220663) q[0];
sx q[0];
rz(2.2882373) q[0];
rz(-2.536643) q[1];
sx q[1];
rz(-1.9974983) q[1];
sx q[1];
rz(-2.6430184) q[1];
rz(2.5685979) q[2];
sx q[2];
rz(-1.6934494) q[2];
sx q[2];
rz(2.1394503) q[2];
rz(2.6400186) q[3];
sx q[3];
rz(-1.4091103) q[3];
sx q[3];
rz(1.4179358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
