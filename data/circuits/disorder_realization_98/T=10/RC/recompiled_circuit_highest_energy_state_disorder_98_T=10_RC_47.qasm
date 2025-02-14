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
rz(0.75818169) q[0];
sx q[0];
rz(3.5092847) q[0];
sx q[0];
rz(7.379471) q[0];
rz(0.81049377) q[1];
sx q[1];
rz(-0.23263045) q[1];
sx q[1];
rz(-2.0356324) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6152412) q[0];
sx q[0];
rz(-0.59068524) q[0];
sx q[0];
rz(0.75449852) q[0];
x q[1];
rz(-1.4742548) q[2];
sx q[2];
rz(-1.8372922) q[2];
sx q[2];
rz(-2.0713553) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.76274058) q[1];
sx q[1];
rz(-1.5334852) q[1];
sx q[1];
rz(0.0040405063) q[1];
x q[2];
rz(0.28438045) q[3];
sx q[3];
rz(-1.5640946) q[3];
sx q[3];
rz(-2.0363765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.37630633) q[2];
sx q[2];
rz(-2.6708965) q[2];
sx q[2];
rz(-2.7619696) q[2];
rz(1.6342573) q[3];
sx q[3];
rz(-1.0999271) q[3];
sx q[3];
rz(-0.96989337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4955502) q[0];
sx q[0];
rz(-2.8046785) q[0];
sx q[0];
rz(-0.90091339) q[0];
rz(0.13149978) q[1];
sx q[1];
rz(-1.200518) q[1];
sx q[1];
rz(0.44201717) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1832804) q[0];
sx q[0];
rz(-1.6526451) q[0];
sx q[0];
rz(0.32572066) q[0];
x q[1];
rz(0.0046185812) q[2];
sx q[2];
rz(-1.5698264) q[2];
sx q[2];
rz(-2.0624954) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9513248) q[1];
sx q[1];
rz(-0.76808483) q[1];
sx q[1];
rz(-0.70835857) q[1];
rz(-pi) q[2];
rz(-1.4591181) q[3];
sx q[3];
rz(-1.1726716) q[3];
sx q[3];
rz(-2.3071837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.76571959) q[2];
sx q[2];
rz(-1.4292382) q[2];
sx q[2];
rz(0.6089375) q[2];
rz(-2.7385312) q[3];
sx q[3];
rz(-1.9654704) q[3];
sx q[3];
rz(-2.6398931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32560638) q[0];
sx q[0];
rz(-2.7222962) q[0];
sx q[0];
rz(1.1380648) q[0];
rz(-1.552938) q[1];
sx q[1];
rz(-0.76858968) q[1];
sx q[1];
rz(0.80702153) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0332542) q[0];
sx q[0];
rz(-1.4716118) q[0];
sx q[0];
rz(1.5434815) q[0];
rz(-2.8093178) q[2];
sx q[2];
rz(-2.6746779) q[2];
sx q[2];
rz(0.98540598) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50259604) q[1];
sx q[1];
rz(-1.6901053) q[1];
sx q[1];
rz(-2.9519777) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.85523) q[3];
sx q[3];
rz(-1.2060646) q[3];
sx q[3];
rz(1.3160365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1300065) q[2];
sx q[2];
rz(-0.61379543) q[2];
sx q[2];
rz(-1.2137132) q[2];
rz(-0.064149292) q[3];
sx q[3];
rz(-1.4870653) q[3];
sx q[3];
rz(0.00042644342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5478058) q[0];
sx q[0];
rz(-1.8812027) q[0];
sx q[0];
rz(2.2547145) q[0];
rz(-2.9828494) q[1];
sx q[1];
rz(-2.6491149) q[1];
sx q[1];
rz(-1.8633206) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87065164) q[0];
sx q[0];
rz(-2.2566911) q[0];
sx q[0];
rz(0.57620462) q[0];
x q[1];
rz(-1.099894) q[2];
sx q[2];
rz(-1.0702025) q[2];
sx q[2];
rz(-2.4125227) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.51047303) q[1];
sx q[1];
rz(-2.2620631) q[1];
sx q[1];
rz(-3.0635628) q[1];
rz(-pi) q[2];
rz(2.9904891) q[3];
sx q[3];
rz(-2.1204815) q[3];
sx q[3];
rz(1.0719887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0025803) q[2];
sx q[2];
rz(-0.85550344) q[2];
sx q[2];
rz(2.1878237) q[2];
rz(-1.194713) q[3];
sx q[3];
rz(-2.298893) q[3];
sx q[3];
rz(0.4755303) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3139528) q[0];
sx q[0];
rz(-2.4251745) q[0];
sx q[0];
rz(-0.59999505) q[0];
rz(-2.4586239) q[1];
sx q[1];
rz(-2.3536847) q[1];
sx q[1];
rz(2.8111828) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67067388) q[0];
sx q[0];
rz(-0.74711159) q[0];
sx q[0];
rz(-1.8404191) q[0];
rz(2.9783713) q[2];
sx q[2];
rz(-1.7931869) q[2];
sx q[2];
rz(1.4537653) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7327795) q[1];
sx q[1];
rz(-2.6169852) q[1];
sx q[1];
rz(1.9636159) q[1];
rz(-1.0903301) q[3];
sx q[3];
rz(-0.76101979) q[3];
sx q[3];
rz(-0.42945592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7052475) q[2];
sx q[2];
rz(-1.0785495) q[2];
sx q[2];
rz(-2.516563) q[2];
rz(2.446512) q[3];
sx q[3];
rz(-1.86097) q[3];
sx q[3];
rz(-0.052791031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.92796749) q[0];
sx q[0];
rz(-1.1518421) q[0];
sx q[0];
rz(1.7684162) q[0];
rz(1.7534509) q[1];
sx q[1];
rz(-2.2699247) q[1];
sx q[1];
rz(1.53481) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2508188) q[0];
sx q[0];
rz(-1.5355331) q[0];
sx q[0];
rz(1.718344) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9589171) q[2];
sx q[2];
rz(-1.2625045) q[2];
sx q[2];
rz(-1.3253554) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.88253792) q[1];
sx q[1];
rz(-1.9758304) q[1];
sx q[1];
rz(-0.71230131) q[1];
rz(-1.9717384) q[3];
sx q[3];
rz(-1.4433268) q[3];
sx q[3];
rz(2.0522224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2024978) q[2];
sx q[2];
rz(-1.9715344) q[2];
sx q[2];
rz(-1.5411752) q[2];
rz(-1.0209068) q[3];
sx q[3];
rz(-1.7424135) q[3];
sx q[3];
rz(-0.30103621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0114667) q[0];
sx q[0];
rz(-0.13700329) q[0];
sx q[0];
rz(2.2504508) q[0];
rz(2.4961684) q[1];
sx q[1];
rz(-1.7452469) q[1];
sx q[1];
rz(2.3088764) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6785936) q[0];
sx q[0];
rz(-0.4834673) q[0];
sx q[0];
rz(0.2966169) q[0];
rz(-pi) q[1];
rz(1.4786167) q[2];
sx q[2];
rz(-0.62609172) q[2];
sx q[2];
rz(0.41958671) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5834076) q[1];
sx q[1];
rz(-1.7654357) q[1];
sx q[1];
rz(-2.6411112) q[1];
rz(-pi) q[2];
rz(-1.1770958) q[3];
sx q[3];
rz(-1.0648921) q[3];
sx q[3];
rz(-0.33123744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5523395) q[2];
sx q[2];
rz(-2.4746042) q[2];
sx q[2];
rz(1.5935295) q[2];
rz(-0.42406905) q[3];
sx q[3];
rz(-1.4196906) q[3];
sx q[3];
rz(0.29485318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.1770723) q[0];
sx q[0];
rz(-1.1966713) q[0];
sx q[0];
rz(0.82343423) q[0];
rz(-1.9442762) q[1];
sx q[1];
rz(-1.6540534) q[1];
sx q[1];
rz(1.6808602) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4890503) q[0];
sx q[0];
rz(-0.33777896) q[0];
sx q[0];
rz(-1.2342288) q[0];
x q[1];
rz(0.97779556) q[2];
sx q[2];
rz(-1.4306746) q[2];
sx q[2];
rz(-1.3205547) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9842661) q[1];
sx q[1];
rz(-1.6132024) q[1];
sx q[1];
rz(2.4799816) q[1];
rz(1.2498871) q[3];
sx q[3];
rz(-1.2276638) q[3];
sx q[3];
rz(0.58716256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0494277) q[2];
sx q[2];
rz(-2.4804513) q[2];
sx q[2];
rz(3.0180422) q[2];
rz(0.83856797) q[3];
sx q[3];
rz(-2.0288012) q[3];
sx q[3];
rz(2.6113094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11414828) q[0];
sx q[0];
rz(-1.0076948) q[0];
sx q[0];
rz(1.4134407) q[0];
rz(2.24276) q[1];
sx q[1];
rz(-2.2347968) q[1];
sx q[1];
rz(-2.5487505) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5699485) q[0];
sx q[0];
rz(-2.1396532) q[0];
sx q[0];
rz(2.4242899) q[0];
rz(-pi) q[1];
rz(2.7024621) q[2];
sx q[2];
rz(-0.20072099) q[2];
sx q[2];
rz(-1.8619271) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.457691) q[1];
sx q[1];
rz(-1.8967034) q[1];
sx q[1];
rz(-2.0052471) q[1];
x q[2];
rz(-2.1778706) q[3];
sx q[3];
rz(-0.98120171) q[3];
sx q[3];
rz(-0.88966767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2841407) q[2];
sx q[2];
rz(-2.393674) q[2];
sx q[2];
rz(0.22135529) q[2];
rz(-2.0981233) q[3];
sx q[3];
rz(-0.9404434) q[3];
sx q[3];
rz(-0.38844696) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8429883) q[0];
sx q[0];
rz(-0.28268155) q[0];
sx q[0];
rz(2.2931732) q[0];
rz(3.0449955) q[1];
sx q[1];
rz(-1.074147) q[1];
sx q[1];
rz(0.75327795) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2601813) q[0];
sx q[0];
rz(-2.8927566) q[0];
sx q[0];
rz(1.9877861) q[0];
rz(-0.033107759) q[2];
sx q[2];
rz(-2.3121142) q[2];
sx q[2];
rz(0.25698369) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.85791432) q[1];
sx q[1];
rz(-1.1925624) q[1];
sx q[1];
rz(2.0973732) q[1];
rz(-2.6274458) q[3];
sx q[3];
rz(-0.78350583) q[3];
sx q[3];
rz(-1.0422106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0775371) q[2];
sx q[2];
rz(-0.83751837) q[2];
sx q[2];
rz(-0.45292863) q[2];
rz(2.5405267) q[3];
sx q[3];
rz(-1.6744924) q[3];
sx q[3];
rz(-0.99630228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30408981) q[0];
sx q[0];
rz(-1.4488198) q[0];
sx q[0];
rz(-2.6813843) q[0];
rz(2.9651463) q[1];
sx q[1];
rz(-2.9063168) q[1];
sx q[1];
rz(1.6520687) q[1];
rz(1.5701736) q[2];
sx q[2];
rz(-1.9194308) q[2];
sx q[2];
rz(0.44375833) q[2];
rz(2.413977) q[3];
sx q[3];
rz(-1.1066827) q[3];
sx q[3];
rz(-0.91409693) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
