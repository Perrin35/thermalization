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
rz(-1.4323956) q[0];
sx q[0];
rz(-0.30269912) q[0];
sx q[0];
rz(-2.1085289) q[0];
rz(1.9167702) q[1];
sx q[1];
rz(1.4900102) q[1];
sx q[1];
rz(9.23041) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8914999) q[0];
sx q[0];
rz(-2.19133) q[0];
sx q[0];
rz(1.3832072) q[0];
rz(1.871045) q[2];
sx q[2];
rz(-0.89648834) q[2];
sx q[2];
rz(-2.6943494) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7558507) q[1];
sx q[1];
rz(-3.0839018) q[1];
sx q[1];
rz(2.2123446) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62567775) q[3];
sx q[3];
rz(-1.2114269) q[3];
sx q[3];
rz(0.43454447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.73244652) q[2];
sx q[2];
rz(-2.0398085) q[2];
sx q[2];
rz(-2.5771602) q[2];
rz(-1.8730646) q[3];
sx q[3];
rz(-3.1334435) q[3];
sx q[3];
rz(0.90659365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.059435189) q[0];
sx q[0];
rz(-3.1340288) q[0];
sx q[0];
rz(-2.2285158) q[0];
rz(0.68436855) q[1];
sx q[1];
rz(-0.00033907779) q[1];
sx q[1];
rz(2.2025462) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8565687) q[0];
sx q[0];
rz(-1.9459581) q[0];
sx q[0];
rz(-0.60174353) q[0];
x q[1];
rz(0.026264391) q[2];
sx q[2];
rz(-0.46411465) q[2];
sx q[2];
rz(0.26192203) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0239183) q[1];
sx q[1];
rz(-0.54421762) q[1];
sx q[1];
rz(-2.2571889) q[1];
rz(-1.7250118) q[3];
sx q[3];
rz(-2.7304478) q[3];
sx q[3];
rz(-1.7257041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3188476) q[2];
sx q[2];
rz(-3.1343967) q[2];
sx q[2];
rz(2.2491271) q[2];
rz(-1.5626296) q[3];
sx q[3];
rz(-3.1171418) q[3];
sx q[3];
rz(-1.310937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7398294) q[0];
sx q[0];
rz(-0.50245291) q[0];
sx q[0];
rz(0.40973642) q[0];
rz(0.01092625) q[1];
sx q[1];
rz(-0.22053638) q[1];
sx q[1];
rz(1.3440557) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7342012) q[0];
sx q[0];
rz(-2.0816861) q[0];
sx q[0];
rz(3.1180918) q[0];
rz(-pi) q[1];
rz(-1.4553237) q[2];
sx q[2];
rz(-1.6342548) q[2];
sx q[2];
rz(0.34380128) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1045055) q[1];
sx q[1];
rz(-1.5059789) q[1];
sx q[1];
rz(-1.3124036) q[1];
rz(1.4859173) q[3];
sx q[3];
rz(-0.1156919) q[3];
sx q[3];
rz(-1.8574238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.672013) q[2];
sx q[2];
rz(-1.6894222) q[2];
sx q[2];
rz(-3.0984042) q[2];
rz(-1.1399266) q[3];
sx q[3];
rz(-2.98525) q[3];
sx q[3];
rz(-1.9388916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-1.7286872) q[0];
sx q[0];
rz(-2.7396956) q[0];
sx q[0];
rz(1.8034978) q[0];
rz(-2.4463704) q[1];
sx q[1];
rz(-3.0137364) q[1];
sx q[1];
rz(-0.41740886) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6175175) q[0];
sx q[0];
rz(-1.39912) q[0];
sx q[0];
rz(2.2508932) q[0];
rz(-pi) q[1];
rz(2.874006) q[2];
sx q[2];
rz(-2.1998465) q[2];
sx q[2];
rz(2.3305394) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7898832) q[1];
sx q[1];
rz(-1.7725866) q[1];
sx q[1];
rz(2.2191597) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8700852) q[3];
sx q[3];
rz(-2.6911754) q[3];
sx q[3];
rz(1.822644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43856105) q[2];
sx q[2];
rz(-2.265354) q[2];
sx q[2];
rz(1.38928) q[2];
rz(-0.82052463) q[3];
sx q[3];
rz(-1.5912143) q[3];
sx q[3];
rz(-2.4337721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.163212) q[0];
sx q[0];
rz(-0.037540171) q[0];
sx q[0];
rz(0.96330825) q[0];
rz(0.18927255) q[1];
sx q[1];
rz(-0.015733868) q[1];
sx q[1];
rz(-0.16545573) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67901299) q[0];
sx q[0];
rz(-1.6143342) q[0];
sx q[0];
rz(1.5971558) q[0];
rz(-2.2654183) q[2];
sx q[2];
rz(-1.3378222) q[2];
sx q[2];
rz(-0.35859749) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9349337) q[1];
sx q[1];
rz(-0.77261415) q[1];
sx q[1];
rz(-1.4040703) q[1];
rz(-1.1375764) q[3];
sx q[3];
rz(-1.8611844) q[3];
sx q[3];
rz(-3.0620787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8339707) q[2];
sx q[2];
rz(-0.51887363) q[2];
sx q[2];
rz(2.0818254) q[2];
rz(-2.4506532) q[3];
sx q[3];
rz(-0.86425224) q[3];
sx q[3];
rz(1.2714269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5237913) q[0];
sx q[0];
rz(-3.0484564) q[0];
sx q[0];
rz(-1.5366489) q[0];
rz(-2.6887584) q[1];
sx q[1];
rz(-0.0080527877) q[1];
sx q[1];
rz(1.4401999) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0217845) q[0];
sx q[0];
rz(-2.9077466) q[0];
sx q[0];
rz(-2.3709557) q[0];
rz(1.790449) q[2];
sx q[2];
rz(-1.7872602) q[2];
sx q[2];
rz(-0.6714657) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3041087) q[1];
sx q[1];
rz(-0.15481259) q[1];
sx q[1];
rz(-0.047708851) q[1];
rz(-pi) q[2];
rz(3.0423445) q[3];
sx q[3];
rz(-1.4008879) q[3];
sx q[3];
rz(2.3836294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3692533) q[2];
sx q[2];
rz(-2.1489096) q[2];
sx q[2];
rz(0.079027979) q[2];
rz(-0.8655656) q[3];
sx q[3];
rz(-2.1871958) q[3];
sx q[3];
rz(-1.0237833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2257776) q[0];
sx q[0];
rz(-3.1362035) q[0];
sx q[0];
rz(2.916577) q[0];
rz(-2.8354538) q[1];
sx q[1];
rz(-3.1248416) q[1];
sx q[1];
rz(-2.2711145) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3358094) q[0];
sx q[0];
rz(-0.085698232) q[0];
sx q[0];
rz(2.5524469) q[0];
rz(-pi) q[1];
rz(-2.4267459) q[2];
sx q[2];
rz(-2.0185592) q[2];
sx q[2];
rz(3.062742) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6806769) q[1];
sx q[1];
rz(-2.3391238) q[1];
sx q[1];
rz(1.6548619) q[1];
rz(-pi) q[2];
rz(-2.2629396) q[3];
sx q[3];
rz(-0.60135287) q[3];
sx q[3];
rz(-0.9200615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.47591448) q[2];
sx q[2];
rz(-0.12683882) q[2];
sx q[2];
rz(1.0352943) q[2];
rz(0.78694844) q[3];
sx q[3];
rz(-3.0988155) q[3];
sx q[3];
rz(2.8662445) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1110474) q[0];
sx q[0];
rz(-0.011140911) q[0];
sx q[0];
rz(3.1159478) q[0];
rz(1.8772839) q[1];
sx q[1];
rz(-0.023921078) q[1];
sx q[1];
rz(-0.69902507) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2866488) q[0];
sx q[0];
rz(-0.35084769) q[0];
sx q[0];
rz(2.4445266) q[0];
x q[1];
rz(1.7412296) q[2];
sx q[2];
rz(-1.1079857) q[2];
sx q[2];
rz(1.5968666) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5311218) q[1];
sx q[1];
rz(-1.4969331) q[1];
sx q[1];
rz(0.38458891) q[1];
x q[2];
rz(2.2243629) q[3];
sx q[3];
rz(-1.1207657) q[3];
sx q[3];
rz(0.91564028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8890624) q[2];
sx q[2];
rz(-2.3928596) q[2];
sx q[2];
rz(-0.32386455) q[2];
rz(-2.9845386) q[3];
sx q[3];
rz(-1.9040949) q[3];
sx q[3];
rz(2.9878555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5667628) q[0];
sx q[0];
rz(-3.1269508) q[0];
sx q[0];
rz(-0.59004849) q[0];
rz(-0.7198965) q[1];
sx q[1];
rz(-3.0820334) q[1];
sx q[1];
rz(0.86729008) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2197207) q[0];
sx q[0];
rz(-2.2014473) q[0];
sx q[0];
rz(1.7008855) q[0];
rz(1.9407104) q[2];
sx q[2];
rz(-0.6354161) q[2];
sx q[2];
rz(-0.70178343) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.54360753) q[1];
sx q[1];
rz(-1.6389009) q[1];
sx q[1];
rz(-1.6833413) q[1];
rz(-0.85764472) q[3];
sx q[3];
rz(-1.4271171) q[3];
sx q[3];
rz(1.194467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1934293) q[2];
sx q[2];
rz(-2.2488504) q[2];
sx q[2];
rz(0.11290045) q[2];
rz(2.0924977) q[3];
sx q[3];
rz(-1.5821404) q[3];
sx q[3];
rz(0.73672867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.076544) q[0];
sx q[0];
rz(-1.5815409) q[0];
sx q[0];
rz(-1.5034058) q[0];
rz(-0.63649559) q[1];
sx q[1];
rz(-2.681585) q[1];
sx q[1];
rz(1.5696625) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29964744) q[0];
sx q[0];
rz(-1.491786) q[0];
sx q[0];
rz(-1.6539198) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1384597) q[2];
sx q[2];
rz(-1.5701446) q[2];
sx q[2];
rz(-2.3359483) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0565785) q[1];
sx q[1];
rz(-0.002925891) q[1];
sx q[1];
rz(0.080119328) q[1];
rz(-pi) q[2];
rz(1.4079553) q[3];
sx q[3];
rz(-2.0800839) q[3];
sx q[3];
rz(-2.7247938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87127176) q[2];
sx q[2];
rz(-1.6297636) q[2];
sx q[2];
rz(-0.16066571) q[2];
rz(-1.9130982) q[3];
sx q[3];
rz(-0.03385032) q[3];
sx q[3];
rz(-2.9401949) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0495618) q[0];
sx q[0];
rz(-1.1783538) q[0];
sx q[0];
rz(0.023068064) q[0];
rz(-1.6231712) q[1];
sx q[1];
rz(-2.7802614) q[1];
sx q[1];
rz(-2.8526715) q[1];
rz(-1.2779909) q[2];
sx q[2];
rz(-1.59272) q[2];
sx q[2];
rz(-1.1772173) q[2];
rz(2.992441) q[3];
sx q[3];
rz(-1.1356529) q[3];
sx q[3];
rz(0.54855772) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
