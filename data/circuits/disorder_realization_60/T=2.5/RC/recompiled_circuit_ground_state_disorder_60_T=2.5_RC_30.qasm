OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.056501) q[0];
sx q[0];
rz(-2.8488475) q[0];
sx q[0];
rz(-2.2226287) q[0];
rz(2.7594944) q[1];
sx q[1];
rz(-0.16799071) q[1];
sx q[1];
rz(-1.1408495) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56132434) q[0];
sx q[0];
rz(-1.7467357) q[0];
sx q[0];
rz(0.12138155) q[0];
rz(-pi) q[1];
rz(1.2873285) q[2];
sx q[2];
rz(-1.2670484) q[2];
sx q[2];
rz(-1.8423353) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6534868) q[1];
sx q[1];
rz(-1.7035023) q[1];
sx q[1];
rz(-0.38487558) q[1];
x q[2];
rz(-0.13762044) q[3];
sx q[3];
rz(-2.0064266) q[3];
sx q[3];
rz(-2.9916258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.69316489) q[2];
sx q[2];
rz(-2.0824671) q[2];
sx q[2];
rz(-1.1761752) q[2];
rz(0.042393427) q[3];
sx q[3];
rz(-1.3375125) q[3];
sx q[3];
rz(0.5347518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.46928826) q[0];
sx q[0];
rz(-2.7870218) q[0];
sx q[0];
rz(2.9571423) q[0];
rz(0.09672673) q[1];
sx q[1];
rz(-1.6177142) q[1];
sx q[1];
rz(-1.2299889) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9365777) q[0];
sx q[0];
rz(-2.7194571) q[0];
sx q[0];
rz(0.5136712) q[0];
x q[1];
rz(-0.076893496) q[2];
sx q[2];
rz(-0.56083365) q[2];
sx q[2];
rz(0.70181134) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0470386) q[1];
sx q[1];
rz(-2.8287005) q[1];
sx q[1];
rz(-2.439586) q[1];
x q[2];
rz(-3.1026479) q[3];
sx q[3];
rz(-1.9488261) q[3];
sx q[3];
rz(1.8824487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.35473287) q[2];
sx q[2];
rz(-2.730098) q[2];
sx q[2];
rz(-3.1301609) q[2];
rz(1.3139906) q[3];
sx q[3];
rz(-1.7766137) q[3];
sx q[3];
rz(0.7029117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3979724) q[0];
sx q[0];
rz(-0.13298661) q[0];
sx q[0];
rz(-2.5308894) q[0];
rz(0.01344219) q[1];
sx q[1];
rz(-1.2862658) q[1];
sx q[1];
rz(-0.59649831) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72815182) q[0];
sx q[0];
rz(-2.1755009) q[0];
sx q[0];
rz(-2.938943) q[0];
rz(-pi) q[1];
rz(1.1223511) q[2];
sx q[2];
rz(-0.62659133) q[2];
sx q[2];
rz(2.5950876) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1787179) q[1];
sx q[1];
rz(-1.8824008) q[1];
sx q[1];
rz(-2.3942663) q[1];
rz(1.0168996) q[3];
sx q[3];
rz(-1.5427329) q[3];
sx q[3];
rz(1.9848167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4639123) q[2];
sx q[2];
rz(-1.8296506) q[2];
sx q[2];
rz(-0.00027351969) q[2];
rz(-1.6655946) q[3];
sx q[3];
rz(-2.4725584) q[3];
sx q[3];
rz(-3.0345548) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45604712) q[0];
sx q[0];
rz(-2.9750329) q[0];
sx q[0];
rz(-1.1628994) q[0];
rz(-0.97745013) q[1];
sx q[1];
rz(-1.5403427) q[1];
sx q[1];
rz(1.5553442) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68846164) q[0];
sx q[0];
rz(-1.5630088) q[0];
sx q[0];
rz(-1.580419) q[0];
rz(-pi) q[1];
rz(-2.3562217) q[2];
sx q[2];
rz(-1.6024688) q[2];
sx q[2];
rz(0.23508628) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0229637) q[1];
sx q[1];
rz(-1.5563772) q[1];
sx q[1];
rz(-2.6447536) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0078251) q[3];
sx q[3];
rz(-2.6053782) q[3];
sx q[3];
rz(1.0325583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.9419452) q[2];
sx q[2];
rz(-0.5117828) q[2];
sx q[2];
rz(0.23435782) q[2];
rz(2.6394081) q[3];
sx q[3];
rz(-2.1000704) q[3];
sx q[3];
rz(2.485399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31768826) q[0];
sx q[0];
rz(-0.81517878) q[0];
sx q[0];
rz(1.3932047) q[0];
rz(-2.4488917) q[1];
sx q[1];
rz(-2.0557978) q[1];
sx q[1];
rz(2.0511973) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4596953) q[0];
sx q[0];
rz(-0.57006627) q[0];
sx q[0];
rz(-1.1718114) q[0];
x q[1];
rz(-0.52509016) q[2];
sx q[2];
rz(-1.0193362) q[2];
sx q[2];
rz(2.9208825) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7140035) q[1];
sx q[1];
rz(-1.98381) q[1];
sx q[1];
rz(2.9772813) q[1];
x q[2];
rz(1.3133009) q[3];
sx q[3];
rz(-2.5867043) q[3];
sx q[3];
rz(-0.99770791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.580487) q[2];
sx q[2];
rz(-1.9876391) q[2];
sx q[2];
rz(0.53606501) q[2];
rz(-2.8481893) q[3];
sx q[3];
rz(-1.9304201) q[3];
sx q[3];
rz(-2.6611879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90477657) q[0];
sx q[0];
rz(-0.95229709) q[0];
sx q[0];
rz(-0.099040898) q[0];
rz(1.0320484) q[1];
sx q[1];
rz(-0.85056225) q[1];
sx q[1];
rz(-1.438407) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28171646) q[0];
sx q[0];
rz(-1.8252354) q[0];
sx q[0];
rz(-2.2734725) q[0];
rz(-1.4969669) q[2];
sx q[2];
rz(-1.8535239) q[2];
sx q[2];
rz(1.7746995) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8107103) q[1];
sx q[1];
rz(-1.5473978) q[1];
sx q[1];
rz(2.3512273) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.766633) q[3];
sx q[3];
rz(-2.8228033) q[3];
sx q[3];
rz(2.6998883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9859887) q[2];
sx q[2];
rz(-2.6376548) q[2];
sx q[2];
rz(2.3206553) q[2];
rz(-1.4486897) q[3];
sx q[3];
rz(-0.71805787) q[3];
sx q[3];
rz(-1.6528355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89707017) q[0];
sx q[0];
rz(-2.2269766) q[0];
sx q[0];
rz(-2.6389417) q[0];
rz(-1.2333168) q[1];
sx q[1];
rz(-0.93524593) q[1];
sx q[1];
rz(0.9800235) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7293046) q[0];
sx q[0];
rz(-1.8192023) q[0];
sx q[0];
rz(-2.8969943) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29666846) q[2];
sx q[2];
rz(-0.7985332) q[2];
sx q[2];
rz(-0.980033) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7240643) q[1];
sx q[1];
rz(-1.8019946) q[1];
sx q[1];
rz(-0.75059143) q[1];
rz(-pi) q[2];
rz(2.3013744) q[3];
sx q[3];
rz(-1.7257236) q[3];
sx q[3];
rz(-1.8904101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9194453) q[2];
sx q[2];
rz(-1.8541226) q[2];
sx q[2];
rz(-1.7670828) q[2];
rz(1.8253271) q[3];
sx q[3];
rz(-0.85597435) q[3];
sx q[3];
rz(0.98178664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.114349) q[0];
sx q[0];
rz(-2.7482432) q[0];
sx q[0];
rz(0.91841665) q[0];
rz(-2.9367327) q[1];
sx q[1];
rz(-1.058895) q[1];
sx q[1];
rz(-1.4835666) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12143887) q[0];
sx q[0];
rz(-0.32628548) q[0];
sx q[0];
rz(1.0980647) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7006364) q[2];
sx q[2];
rz(-2.3774638) q[2];
sx q[2];
rz(-1.4997327) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0369954) q[1];
sx q[1];
rz(-0.87338398) q[1];
sx q[1];
rz(-3.0265977) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5213883) q[3];
sx q[3];
rz(-0.8522343) q[3];
sx q[3];
rz(-0.75203943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.001699) q[2];
sx q[2];
rz(-0.49109444) q[2];
sx q[2];
rz(-1.8074544) q[2];
rz(0.88207465) q[3];
sx q[3];
rz(-0.69552723) q[3];
sx q[3];
rz(1.849256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(3.0598711) q[0];
sx q[0];
rz(-0.43161714) q[0];
sx q[0];
rz(0.01817848) q[0];
rz(-2.7320618) q[1];
sx q[1];
rz(-1.7117701) q[1];
sx q[1];
rz(3.0407564) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6438599) q[0];
sx q[0];
rz(-0.63526216) q[0];
sx q[0];
rz(-0.058202581) q[0];
rz(-pi) q[1];
rz(1.5201496) q[2];
sx q[2];
rz(-1.1505923) q[2];
sx q[2];
rz(-0.24383185) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3506895) q[1];
sx q[1];
rz(-1.8115582) q[1];
sx q[1];
rz(1.945182) q[1];
rz(1.1222029) q[3];
sx q[3];
rz(-0.27315419) q[3];
sx q[3];
rz(-2.7970683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.33710256) q[2];
sx q[2];
rz(-3.0323995) q[2];
sx q[2];
rz(-1.2479372) q[2];
rz(1.9167871) q[3];
sx q[3];
rz(-0.8453415) q[3];
sx q[3];
rz(3.0758408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(2.9058022) q[0];
sx q[0];
rz(-1.6658655) q[0];
sx q[0];
rz(2.8210848) q[0];
rz(2.7505752) q[1];
sx q[1];
rz(-2.0000439) q[1];
sx q[1];
rz(-1.5919707) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3452197) q[0];
sx q[0];
rz(-1.380305) q[0];
sx q[0];
rz(2.9389771) q[0];
x q[1];
rz(-1.7334601) q[2];
sx q[2];
rz(-1.062359) q[2];
sx q[2];
rz(2.4738653) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.69471162) q[1];
sx q[1];
rz(-1.9384192) q[1];
sx q[1];
rz(3.047961) q[1];
rz(-pi) q[2];
rz(2.7192468) q[3];
sx q[3];
rz(-1.5776488) q[3];
sx q[3];
rz(2.6391878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0512507) q[2];
sx q[2];
rz(-1.0498468) q[2];
sx q[2];
rz(-0.40038294) q[2];
rz(-2.9054437) q[3];
sx q[3];
rz(-0.103424) q[3];
sx q[3];
rz(-1.0108488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(0.26265963) q[0];
sx q[0];
rz(-2.6317609) q[0];
sx q[0];
rz(1.2608933) q[0];
rz(-0.48514584) q[1];
sx q[1];
rz(-0.79882516) q[1];
sx q[1];
rz(-0.47732236) q[1];
rz(1.0322124) q[2];
sx q[2];
rz(-0.24145024) q[2];
sx q[2];
rz(-0.35035124) q[2];
rz(1.3427721) q[3];
sx q[3];
rz(-1.6327471) q[3];
sx q[3];
rz(-3.0625797) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
