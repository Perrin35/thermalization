OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.085091703) q[0];
sx q[0];
rz(-0.29274517) q[0];
sx q[0];
rz(2.2226287) q[0];
rz(-0.38209823) q[1];
sx q[1];
rz(-2.9736019) q[1];
sx q[1];
rz(1.1408495) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5802683) q[0];
sx q[0];
rz(-1.394857) q[0];
sx q[0];
rz(-0.12138155) q[0];
rz(-pi) q[1];
rz(0.31556712) q[2];
sx q[2];
rz(-1.3006388) q[2];
sx q[2];
rz(-0.18462727) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1124494) q[1];
sx q[1];
rz(-1.9521128) q[1];
sx q[1];
rz(-1.7138395) q[1];
rz(-pi) q[2];
rz(-1.2841746) q[3];
sx q[3];
rz(-0.45551963) q[3];
sx q[3];
rz(-2.974433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.69316489) q[2];
sx q[2];
rz(-1.0591256) q[2];
sx q[2];
rz(-1.9654174) q[2];
rz(3.0991992) q[3];
sx q[3];
rz(-1.3375125) q[3];
sx q[3];
rz(-0.5347518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.6723044) q[0];
sx q[0];
rz(-2.7870218) q[0];
sx q[0];
rz(0.18445036) q[0];
rz(-3.0448659) q[1];
sx q[1];
rz(-1.5238785) q[1];
sx q[1];
rz(1.2299889) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.205015) q[0];
sx q[0];
rz(-0.42213556) q[0];
sx q[0];
rz(2.6279215) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.582091) q[2];
sx q[2];
rz(-1.6116665) q[2];
sx q[2];
rz(-0.80383435) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1538196) q[1];
sx q[1];
rz(-1.3706939) q[1];
sx q[1];
rz(0.2421744) q[1];
x q[2];
rz(-1.1925063) q[3];
sx q[3];
rz(-1.6069901) q[3];
sx q[3];
rz(-0.32603273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7868598) q[2];
sx q[2];
rz(-0.41149461) q[2];
sx q[2];
rz(-3.1301609) q[2];
rz(1.8276021) q[3];
sx q[3];
rz(-1.3649789) q[3];
sx q[3];
rz(0.7029117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74362022) q[0];
sx q[0];
rz(-3.008606) q[0];
sx q[0];
rz(-2.5308894) q[0];
rz(-0.01344219) q[1];
sx q[1];
rz(-1.8553269) q[1];
sx q[1];
rz(2.5450943) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0666351) q[0];
sx q[0];
rz(-0.63371113) q[0];
sx q[0];
rz(1.8541965) q[0];
x q[1];
rz(0.30412401) q[2];
sx q[2];
rz(-2.1275188) q[2];
sx q[2];
rz(-3.1310658) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1787179) q[1];
sx q[1];
rz(-1.8824008) q[1];
sx q[1];
rz(-2.3942663) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0168996) q[3];
sx q[3];
rz(-1.5427329) q[3];
sx q[3];
rz(-1.9848167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.67768031) q[2];
sx q[2];
rz(-1.8296506) q[2];
sx q[2];
rz(0.00027351969) q[2];
rz(1.6655946) q[3];
sx q[3];
rz(-2.4725584) q[3];
sx q[3];
rz(3.0345548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6855455) q[0];
sx q[0];
rz(-0.16655971) q[0];
sx q[0];
rz(-1.1628994) q[0];
rz(-0.97745013) q[1];
sx q[1];
rz(-1.5403427) q[1];
sx q[1];
rz(1.5553442) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88240963) q[0];
sx q[0];
rz(-1.5804187) q[0];
sx q[0];
rz(-0.0077879328) q[0];
rz(-pi) q[1];
rz(2.3562217) q[2];
sx q[2];
rz(-1.5391239) q[2];
sx q[2];
rz(-2.9065064) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4443497) q[1];
sx q[1];
rz(-2.067579) q[1];
sx q[1];
rz(1.5871983) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0647207) q[3];
sx q[3];
rz(-1.7887497) q[3];
sx q[3];
rz(0.92011425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1996475) q[2];
sx q[2];
rz(-2.6298099) q[2];
sx q[2];
rz(0.23435782) q[2];
rz(0.50218454) q[3];
sx q[3];
rz(-2.1000704) q[3];
sx q[3];
rz(-2.485399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8239044) q[0];
sx q[0];
rz(-0.81517878) q[0];
sx q[0];
rz(-1.748388) q[0];
rz(-0.69270095) q[1];
sx q[1];
rz(-2.0557978) q[1];
sx q[1];
rz(1.0903953) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2176184) q[0];
sx q[0];
rz(-1.0503142) q[0];
sx q[0];
rz(0.24407669) q[0];
rz(2.6165025) q[2];
sx q[2];
rz(-1.0193362) q[2];
sx q[2];
rz(2.9208825) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.322256) q[1];
sx q[1];
rz(-2.6988479) q[1];
sx q[1];
rz(-1.2135452) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8282918) q[3];
sx q[3];
rz(-2.5867043) q[3];
sx q[3];
rz(2.1438847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.56110567) q[2];
sx q[2];
rz(-1.1539536) q[2];
sx q[2];
rz(-0.53606501) q[2];
rz(0.29340336) q[3];
sx q[3];
rz(-1.9304201) q[3];
sx q[3];
rz(-2.6611879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90477657) q[0];
sx q[0];
rz(-2.1892956) q[0];
sx q[0];
rz(-0.099040898) q[0];
rz(-2.1095443) q[1];
sx q[1];
rz(-2.2910304) q[1];
sx q[1];
rz(-1.7031857) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28171646) q[0];
sx q[0];
rz(-1.8252354) q[0];
sx q[0];
rz(2.2734725) q[0];
x q[1];
rz(2.8929404) q[2];
sx q[2];
rz(-0.2919582) q[2];
sx q[2];
rz(-2.0338634) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26353961) q[1];
sx q[1];
rz(-0.78070736) q[1];
sx q[1];
rz(1.5375464) q[1];
rz(-1.4505054) q[3];
sx q[3];
rz(-1.2748537) q[3];
sx q[3];
rz(3.0927998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9859887) q[2];
sx q[2];
rz(-2.6376548) q[2];
sx q[2];
rz(2.3206553) q[2];
rz(1.692903) q[3];
sx q[3];
rz(-2.4235348) q[3];
sx q[3];
rz(-1.4887571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89707017) q[0];
sx q[0];
rz(-2.2269766) q[0];
sx q[0];
rz(-0.50265092) q[0];
rz(1.9082759) q[1];
sx q[1];
rz(-0.93524593) q[1];
sx q[1];
rz(-2.1615692) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4122881) q[0];
sx q[0];
rz(-1.8192023) q[0];
sx q[0];
rz(-0.24459837) q[0];
rz(0.77620164) q[2];
sx q[2];
rz(-1.3598249) q[2];
sx q[2];
rz(0.38061505) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.087638559) q[1];
sx q[1];
rz(-0.77869895) q[1];
sx q[1];
rz(2.8092572) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9348922) q[3];
sx q[3];
rz(-2.2906942) q[3];
sx q[3];
rz(-2.9593619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.22214733) q[2];
sx q[2];
rz(-1.8541226) q[2];
sx q[2];
rz(1.3745098) q[2];
rz(-1.8253271) q[3];
sx q[3];
rz(-0.85597435) q[3];
sx q[3];
rz(-0.98178664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027243622) q[0];
sx q[0];
rz(-2.7482432) q[0];
sx q[0];
rz(-2.223176) q[0];
rz(2.9367327) q[1];
sx q[1];
rz(-2.0826976) q[1];
sx q[1];
rz(1.6580261) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0201538) q[0];
sx q[0];
rz(-0.32628548) q[0];
sx q[0];
rz(1.0980647) q[0];
rz(2.3307033) q[2];
sx q[2];
rz(-1.6605018) q[2];
sx q[2];
rz(-2.9765338) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.073347884) q[1];
sx q[1];
rz(-2.4363344) q[1];
sx q[1];
rz(1.4346992) q[1];
rz(-pi) q[2];
rz(3.0851787) q[3];
sx q[3];
rz(-0.71995633) q[3];
sx q[3];
rz(2.3145793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.001699) q[2];
sx q[2];
rz(-0.49109444) q[2];
sx q[2];
rz(1.8074544) q[2];
rz(0.88207465) q[3];
sx q[3];
rz(-2.4460654) q[3];
sx q[3];
rz(-1.849256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081721574) q[0];
sx q[0];
rz(-2.7099755) q[0];
sx q[0];
rz(3.1234142) q[0];
rz(-0.40953088) q[1];
sx q[1];
rz(-1.7117701) q[1];
sx q[1];
rz(0.10083625) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6438599) q[0];
sx q[0];
rz(-2.5063305) q[0];
sx q[0];
rz(-3.0833901) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5201496) q[2];
sx q[2];
rz(-1.1505923) q[2];
sx q[2];
rz(2.8977608) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8162569) q[1];
sx q[1];
rz(-0.44202572) q[1];
sx q[1];
rz(0.97953743) q[1];
rz(-1.8180679) q[3];
sx q[3];
rz(-1.4535289) q[3];
sx q[3];
rz(2.3493249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.33710256) q[2];
sx q[2];
rz(-3.0323995) q[2];
sx q[2];
rz(1.8936554) q[2];
rz(1.9167871) q[3];
sx q[3];
rz(-2.2962511) q[3];
sx q[3];
rz(-3.0758408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2357904) q[0];
sx q[0];
rz(-1.4757272) q[0];
sx q[0];
rz(-2.8210848) q[0];
rz(-2.7505752) q[1];
sx q[1];
rz(-1.1415488) q[1];
sx q[1];
rz(-1.5919707) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3452197) q[0];
sx q[0];
rz(-1.7612877) q[0];
sx q[0];
rz(2.9389771) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28279998) q[2];
sx q[2];
rz(-0.53164266) q[2];
sx q[2];
rz(-0.34257364) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.446881) q[1];
sx q[1];
rz(-1.2031735) q[1];
sx q[1];
rz(0.093631677) q[1];
x q[2];
rz(1.5783089) q[3];
sx q[3];
rz(-1.9931317) q[3];
sx q[3];
rz(-2.0762805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0512507) q[2];
sx q[2];
rz(-1.0498468) q[2];
sx q[2];
rz(-2.7412097) q[2];
rz(-2.9054437) q[3];
sx q[3];
rz(-0.103424) q[3];
sx q[3];
rz(-1.0108488) q[3];
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
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.878933) q[0];
sx q[0];
rz(-2.6317609) q[0];
sx q[0];
rz(1.2608933) q[0];
rz(-0.48514584) q[1];
sx q[1];
rz(-0.79882516) q[1];
sx q[1];
rz(-0.47732236) q[1];
rz(-2.1093802) q[2];
sx q[2];
rz(-0.24145024) q[2];
sx q[2];
rz(-0.35035124) q[2];
rz(-1.3029836) q[3];
sx q[3];
rz(-0.23614863) q[3];
sx q[3];
rz(1.9105259) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
