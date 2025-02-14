OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0786809) q[0];
sx q[0];
rz(-0.39830387) q[0];
sx q[0];
rz(2.8226573) q[0];
rz(-0.4660663) q[1];
sx q[1];
rz(-1.708344) q[1];
sx q[1];
rz(0.27369174) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6382054) q[0];
sx q[0];
rz(-0.84716684) q[0];
sx q[0];
rz(2.6658325) q[0];
rz(2.3897522) q[2];
sx q[2];
rz(-2.813857) q[2];
sx q[2];
rz(0.40204266) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.92299709) q[1];
sx q[1];
rz(-1.6875008) q[1];
sx q[1];
rz(1.7497108) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6998547) q[3];
sx q[3];
rz(-2.3328247) q[3];
sx q[3];
rz(-1.3490775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29229257) q[2];
sx q[2];
rz(-2.2165074) q[2];
sx q[2];
rz(-1.7577457) q[2];
rz(-0.84550953) q[3];
sx q[3];
rz(-0.011761646) q[3];
sx q[3];
rz(-2.1498634) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94075769) q[0];
sx q[0];
rz(-2.2469914) q[0];
sx q[0];
rz(1.0351329) q[0];
rz(-1.9738522) q[1];
sx q[1];
rz(-0.20226856) q[1];
sx q[1];
rz(2.367173) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7073785) q[0];
sx q[0];
rz(-1.0814572) q[0];
sx q[0];
rz(-2.6721832) q[0];
x q[1];
rz(0.34458745) q[2];
sx q[2];
rz(-0.93898749) q[2];
sx q[2];
rz(2.5528204) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.064146994) q[1];
sx q[1];
rz(-1.059491) q[1];
sx q[1];
rz(3.0102986) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6210591) q[3];
sx q[3];
rz(-1.5368965) q[3];
sx q[3];
rz(-0.14824319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3356129) q[2];
sx q[2];
rz(-0.077294417) q[2];
sx q[2];
rz(0.94266164) q[2];
rz(-0.089426905) q[3];
sx q[3];
rz(-2.4783897) q[3];
sx q[3];
rz(0.30557835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-0.59219229) q[0];
sx q[0];
rz(-2.5553199) q[0];
sx q[0];
rz(-0.19202448) q[0];
rz(2.5281455) q[1];
sx q[1];
rz(-2.6934721) q[1];
sx q[1];
rz(-3.1268069) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7887869) q[0];
sx q[0];
rz(-2.3067622) q[0];
sx q[0];
rz(-1.0218746) q[0];
x q[1];
rz(-2.4528265) q[2];
sx q[2];
rz(-1.7851794) q[2];
sx q[2];
rz(-0.87300473) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.80577981) q[1];
sx q[1];
rz(-1.8034013) q[1];
sx q[1];
rz(-0.97877731) q[1];
rz(-pi) q[2];
rz(-1.5824541) q[3];
sx q[3];
rz(-1.5069256) q[3];
sx q[3];
rz(-3.0272878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.52745596) q[2];
sx q[2];
rz(-1.5469896) q[2];
sx q[2];
rz(3.1318437) q[2];
rz(2.5629432) q[3];
sx q[3];
rz(-2.0895683) q[3];
sx q[3];
rz(-0.78154045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11827271) q[0];
sx q[0];
rz(-0.8096205) q[0];
sx q[0];
rz(-0.6672346) q[0];
rz(1.0046593) q[1];
sx q[1];
rz(-2.9861082) q[1];
sx q[1];
rz(-1.3803233) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28128703) q[0];
sx q[0];
rz(-1.5800887) q[0];
sx q[0];
rz(-2.0486661) q[0];
x q[1];
rz(-2.818872) q[2];
sx q[2];
rz(-1.1185794) q[2];
sx q[2];
rz(1.9911901) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2365773) q[1];
sx q[1];
rz(-1.5386174) q[1];
sx q[1];
rz(-2.6479812) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9033086) q[3];
sx q[3];
rz(-1.0589979) q[3];
sx q[3];
rz(-0.82945591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.37463793) q[2];
sx q[2];
rz(-1.1344323) q[2];
sx q[2];
rz(3.0881506) q[2];
rz(2.5080569) q[3];
sx q[3];
rz(-2.7233248) q[3];
sx q[3];
rz(0.0088648349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87090129) q[0];
sx q[0];
rz(-0.57796657) q[0];
sx q[0];
rz(2.0722678) q[0];
rz(1.2879734) q[1];
sx q[1];
rz(-0.063871495) q[1];
sx q[1];
rz(-3.0499444) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1250138) q[0];
sx q[0];
rz(-1.5591987) q[0];
sx q[0];
rz(1.5484823) q[0];
x q[1];
rz(-2.0266125) q[2];
sx q[2];
rz(-1.2741003) q[2];
sx q[2];
rz(1.6216344) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5875467) q[1];
sx q[1];
rz(-0.57922208) q[1];
sx q[1];
rz(-0.86279561) q[1];
rz(-pi) q[2];
rz(2.3717068) q[3];
sx q[3];
rz(-1.5987694) q[3];
sx q[3];
rz(-0.5045435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2962467) q[2];
sx q[2];
rz(-2.3643934) q[2];
sx q[2];
rz(0.11543154) q[2];
rz(0.7655862) q[3];
sx q[3];
rz(-1.4976394) q[3];
sx q[3];
rz(-2.7287741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5033443) q[0];
sx q[0];
rz(-2.4920576) q[0];
sx q[0];
rz(0.58393884) q[0];
rz(-3.0931547) q[1];
sx q[1];
rz(-0.2225288) q[1];
sx q[1];
rz(0.64360523) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22065255) q[0];
sx q[0];
rz(-0.31198129) q[0];
sx q[0];
rz(1.4156154) q[0];
rz(-pi) q[1];
rz(-2.2320643) q[2];
sx q[2];
rz(-1.1937181) q[2];
sx q[2];
rz(2.909735) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.076526) q[1];
sx q[1];
rz(-2.4187633) q[1];
sx q[1];
rz(1.3704872) q[1];
rz(-pi) q[2];
rz(0.77966431) q[3];
sx q[3];
rz(-1.7105967) q[3];
sx q[3];
rz(1.1598905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.863997) q[2];
sx q[2];
rz(-2.3433351) q[2];
sx q[2];
rz(2.4994728) q[2];
rz(0.13907214) q[3];
sx q[3];
rz(-3.0032447) q[3];
sx q[3];
rz(1.4596435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7794401) q[0];
sx q[0];
rz(-2.7208949) q[0];
sx q[0];
rz(0.62533373) q[0];
rz(-0.23896898) q[1];
sx q[1];
rz(-0.23663722) q[1];
sx q[1];
rz(1.3561358) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55955333) q[0];
sx q[0];
rz(-1.6252717) q[0];
sx q[0];
rz(0.031886851) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4932213) q[2];
sx q[2];
rz(-0.4604606) q[2];
sx q[2];
rz(-2.5878276) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0812057) q[1];
sx q[1];
rz(-1.9293509) q[1];
sx q[1];
rz(-0.21292673) q[1];
rz(-2.5111968) q[3];
sx q[3];
rz(-0.26132628) q[3];
sx q[3];
rz(3.0898409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.29505342) q[2];
sx q[2];
rz(-1.2650547) q[2];
sx q[2];
rz(-0.96578252) q[2];
rz(0.24671181) q[3];
sx q[3];
rz(-1.2652946) q[3];
sx q[3];
rz(-1.770796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63916373) q[0];
sx q[0];
rz(-2.736709) q[0];
sx q[0];
rz(0.13614458) q[0];
rz(2.4038521) q[1];
sx q[1];
rz(-2.9164011) q[1];
sx q[1];
rz(1.9053316) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8566187) q[0];
sx q[0];
rz(-2.1559733) q[0];
sx q[0];
rz(-0.43584688) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9604016) q[2];
sx q[2];
rz(-0.94215067) q[2];
sx q[2];
rz(1.9235856) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.63041362) q[1];
sx q[1];
rz(-1.2949589) q[1];
sx q[1];
rz(-2.7008369) q[1];
x q[2];
rz(-2.8131139) q[3];
sx q[3];
rz(-0.94624472) q[3];
sx q[3];
rz(0.49027944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4514734) q[2];
sx q[2];
rz(-1.5855007) q[2];
sx q[2];
rz(-0.95009178) q[2];
rz(-0.39725605) q[3];
sx q[3];
rz(-0.49644956) q[3];
sx q[3];
rz(-0.44601405) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5815247) q[0];
sx q[0];
rz(-0.6811322) q[0];
sx q[0];
rz(-3.0233622) q[0];
rz(1.0826348) q[1];
sx q[1];
rz(-1.2018459) q[1];
sx q[1];
rz(-1.3523678) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23516857) q[0];
sx q[0];
rz(-1.4682859) q[0];
sx q[0];
rz(-0.39484039) q[0];
rz(2.9832411) q[2];
sx q[2];
rz(-2.5917946) q[2];
sx q[2];
rz(-2.9937772) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.39555672) q[1];
sx q[1];
rz(-2.6415351) q[1];
sx q[1];
rz(-1.8386503) q[1];
x q[2];
rz(2.1264137) q[3];
sx q[3];
rz(-2.5930282) q[3];
sx q[3];
rz(-3.0231904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3236986) q[2];
sx q[2];
rz(-2.4511621) q[2];
sx q[2];
rz(-0.77198088) q[2];
rz(0.082993232) q[3];
sx q[3];
rz(-1.8738184) q[3];
sx q[3];
rz(-2.6193589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25987396) q[0];
sx q[0];
rz(-2.3840955) q[0];
sx q[0];
rz(-0.74257332) q[0];
rz(1.8450129) q[1];
sx q[1];
rz(-0.80755889) q[1];
sx q[1];
rz(0.47320941) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5294327) q[0];
sx q[0];
rz(-1.4845779) q[0];
sx q[0];
rz(1.3032806) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80188216) q[2];
sx q[2];
rz(-0.3496162) q[2];
sx q[2];
rz(-2.5494818) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0124346) q[1];
sx q[1];
rz(-1.9170463) q[1];
sx q[1];
rz(1.0770256) q[1];
rz(-pi) q[2];
rz(-2.8678611) q[3];
sx q[3];
rz(-2.377284) q[3];
sx q[3];
rz(-2.4814062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5903198) q[2];
sx q[2];
rz(-1.7251622) q[2];
sx q[2];
rz(1.9080706) q[2];
rz(2.7963855) q[3];
sx q[3];
rz(-2.749741) q[3];
sx q[3];
rz(-0.83693081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45967669) q[0];
sx q[0];
rz(-1.3688594) q[0];
sx q[0];
rz(-1.3874227) q[0];
rz(-0.058902901) q[1];
sx q[1];
rz(-1.420493) q[1];
sx q[1];
rz(-2.4891985) q[1];
rz(1.1299533) q[2];
sx q[2];
rz(-1.9222353) q[2];
sx q[2];
rz(2.0735379) q[2];
rz(0.4555489) q[3];
sx q[3];
rz(-2.1524317) q[3];
sx q[3];
rz(0.47453415) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
