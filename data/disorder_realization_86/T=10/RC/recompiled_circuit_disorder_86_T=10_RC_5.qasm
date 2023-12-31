OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3671626) q[0];
sx q[0];
rz(-2.2280333) q[0];
sx q[0];
rz(1.7295184) q[0];
rz(-2.9867759) q[1];
sx q[1];
rz(5.6875416) q[1];
sx q[1];
rz(7.7653801) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5858551) q[0];
sx q[0];
rz(-1.9940388) q[0];
sx q[0];
rz(-1.4002355) q[0];
rz(-pi) q[1];
rz(2.6719195) q[2];
sx q[2];
rz(-0.28684068) q[2];
sx q[2];
rz(1.6490205) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9285779) q[1];
sx q[1];
rz(-1.9933812) q[1];
sx q[1];
rz(1.0282474) q[1];
rz(3.0406038) q[3];
sx q[3];
rz(-1.0102934) q[3];
sx q[3];
rz(-1.3274173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1564864) q[2];
sx q[2];
rz(-0.50922314) q[2];
sx q[2];
rz(0.86581725) q[2];
rz(-2.1872897) q[3];
sx q[3];
rz(-1.6031957) q[3];
sx q[3];
rz(1.2877864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1433379) q[0];
sx q[0];
rz(-1.4366432) q[0];
sx q[0];
rz(-3.1153733) q[0];
rz(1.6014618) q[1];
sx q[1];
rz(-1.5988348) q[1];
sx q[1];
rz(-2.1781133) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026982633) q[0];
sx q[0];
rz(-0.61404213) q[0];
sx q[0];
rz(1.5747889) q[0];
rz(-2.8112667) q[2];
sx q[2];
rz(-2.1147554) q[2];
sx q[2];
rz(-0.75737539) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0986833) q[1];
sx q[1];
rz(-2.5334362) q[1];
sx q[1];
rz(-0.98867464) q[1];
x q[2];
rz(1.9217334) q[3];
sx q[3];
rz(-1.3388472) q[3];
sx q[3];
rz(-1.9333145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5144689) q[2];
sx q[2];
rz(-1.1273948) q[2];
sx q[2];
rz(-0.13452402) q[2];
rz(-0.7450122) q[3];
sx q[3];
rz(-0.22694215) q[3];
sx q[3];
rz(-0.9427332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9298252) q[0];
sx q[0];
rz(-2.7524502) q[0];
sx q[0];
rz(-0.79743687) q[0];
rz(2.0939317) q[1];
sx q[1];
rz(-0.14973775) q[1];
sx q[1];
rz(2.581596) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13941923) q[0];
sx q[0];
rz(-1.1886485) q[0];
sx q[0];
rz(-0.21811534) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30324869) q[2];
sx q[2];
rz(-1.5504642) q[2];
sx q[2];
rz(3.0603527) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3126038) q[1];
sx q[1];
rz(-1.3723515) q[1];
sx q[1];
rz(0.36638422) q[1];
x q[2];
rz(2.6651354) q[3];
sx q[3];
rz(-1.1343079) q[3];
sx q[3];
rz(-2.1863294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3893163) q[2];
sx q[2];
rz(-1.9439149) q[2];
sx q[2];
rz(-2.9690572) q[2];
rz(-0.98207384) q[3];
sx q[3];
rz(-1.7445824) q[3];
sx q[3];
rz(1.0579695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.003222) q[0];
sx q[0];
rz(-2.0439742) q[0];
sx q[0];
rz(-2.8570783) q[0];
rz(2.8248887) q[1];
sx q[1];
rz(-2.7088294) q[1];
sx q[1];
rz(1.2987312) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5772229) q[0];
sx q[0];
rz(-1.7958613) q[0];
sx q[0];
rz(-1.0610915) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80438517) q[2];
sx q[2];
rz(-1.5717236) q[2];
sx q[2];
rz(1.550012) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5609834) q[1];
sx q[1];
rz(-1.0099851) q[1];
sx q[1];
rz(2.9210655) q[1];
rz(-0.11233791) q[3];
sx q[3];
rz(-1.2445407) q[3];
sx q[3];
rz(2.5374967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6056885) q[2];
sx q[2];
rz(-2.9457592) q[2];
sx q[2];
rz(2.7569125) q[2];
rz(-0.7540594) q[3];
sx q[3];
rz(-1.056517) q[3];
sx q[3];
rz(-1.6872905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6032747) q[0];
sx q[0];
rz(-2.1544927) q[0];
sx q[0];
rz(-1.7549365) q[0];
rz(2.9105913) q[1];
sx q[1];
rz(-1.8004386) q[1];
sx q[1];
rz(0.2968266) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7291527) q[0];
sx q[0];
rz(-1.531633) q[0];
sx q[0];
rz(-2.5705283) q[0];
x q[1];
rz(-2.3987531) q[2];
sx q[2];
rz(-1.4577216) q[2];
sx q[2];
rz(2.6807221) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0823114) q[1];
sx q[1];
rz(-1.1059522) q[1];
sx q[1];
rz(0.44537284) q[1];
rz(-pi) q[2];
rz(0.71367587) q[3];
sx q[3];
rz(-1.6380966) q[3];
sx q[3];
rz(2.0590559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.37830535) q[2];
sx q[2];
rz(-1.8313235) q[2];
sx q[2];
rz(2.7491167) q[2];
rz(1.1522419) q[3];
sx q[3];
rz(-0.71458721) q[3];
sx q[3];
rz(0.31744441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1257989) q[0];
sx q[0];
rz(-1.5725461) q[0];
sx q[0];
rz(0.75138599) q[0];
rz(1.3279351) q[1];
sx q[1];
rz(-1.2633879) q[1];
sx q[1];
rz(0.60633916) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4088926) q[0];
sx q[0];
rz(-3.0928287) q[0];
sx q[0];
rz(2.7932037) q[0];
rz(-pi) q[1];
rz(2.0727856) q[2];
sx q[2];
rz(-0.74337465) q[2];
sx q[2];
rz(0.547264) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.11583466) q[1];
sx q[1];
rz(-2.5587213) q[1];
sx q[1];
rz(1.9028266) q[1];
x q[2];
rz(1.8984406) q[3];
sx q[3];
rz(-2.9122105) q[3];
sx q[3];
rz(-2.4250507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5027344) q[2];
sx q[2];
rz(-1.0417465) q[2];
sx q[2];
rz(1.139337) q[2];
rz(-1.4849439) q[3];
sx q[3];
rz(-1.9610201) q[3];
sx q[3];
rz(3.0373354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5320324) q[0];
sx q[0];
rz(-2.39344) q[0];
sx q[0];
rz(0.50810057) q[0];
rz(-1.5628901) q[1];
sx q[1];
rz(-2.0527614) q[1];
sx q[1];
rz(2.3513444) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7471874) q[0];
sx q[0];
rz(-0.98441511) q[0];
sx q[0];
rz(1.2140843) q[0];
rz(-pi) q[1];
rz(2.9279518) q[2];
sx q[2];
rz(-1.5593312) q[2];
sx q[2];
rz(1.2917047) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7888469) q[1];
sx q[1];
rz(-0.48301304) q[1];
sx q[1];
rz(0.61139815) q[1];
x q[2];
rz(1.1931476) q[3];
sx q[3];
rz(-2.1834063) q[3];
sx q[3];
rz(0.8876422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.053085176) q[2];
sx q[2];
rz(-0.44184703) q[2];
sx q[2];
rz(1.7283758) q[2];
rz(1.6648071) q[3];
sx q[3];
rz(-2.1063185) q[3];
sx q[3];
rz(-3.0800381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7168032) q[0];
sx q[0];
rz(-3.1112818) q[0];
sx q[0];
rz(-2.0943663) q[0];
rz(-0.60910243) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(-1.3887127) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7097276) q[0];
sx q[0];
rz(-2.1301529) q[0];
sx q[0];
rz(3.1064242) q[0];
rz(2.4829364) q[2];
sx q[2];
rz(-2.5052862) q[2];
sx q[2];
rz(3.0898526) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.74275201) q[1];
sx q[1];
rz(-1.5594149) q[1];
sx q[1];
rz(0.091817261) q[1];
rz(-2.8453313) q[3];
sx q[3];
rz(-0.98516609) q[3];
sx q[3];
rz(1.1405917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9528815) q[2];
sx q[2];
rz(-0.41023508) q[2];
sx q[2];
rz(2.2593373) q[2];
rz(-1.7404209) q[3];
sx q[3];
rz(-1.975235) q[3];
sx q[3];
rz(-1.9410979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27154487) q[0];
sx q[0];
rz(-0.4168059) q[0];
sx q[0];
rz(1.7154988) q[0];
rz(3.0601314) q[1];
sx q[1];
rz(-1.1625682) q[1];
sx q[1];
rz(2.5833599) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0575858) q[0];
sx q[0];
rz(-1.1450197) q[0];
sx q[0];
rz(-0.42752479) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.105537) q[2];
sx q[2];
rz(-3.0214587) q[2];
sx q[2];
rz(2.5533954) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3840752) q[1];
sx q[1];
rz(-1.5069403) q[1];
sx q[1];
rz(-2.3149895) q[1];
rz(-1.7562808) q[3];
sx q[3];
rz(-1.0914601) q[3];
sx q[3];
rz(-2.6244147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8490863) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(2.9294087) q[2];
rz(0.21197453) q[3];
sx q[3];
rz(-2.4583355) q[3];
sx q[3];
rz(1.9395444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6367209) q[0];
sx q[0];
rz(-2.3175406) q[0];
sx q[0];
rz(-1.5378392) q[0];
rz(-0.82540712) q[1];
sx q[1];
rz(-0.67276612) q[1];
sx q[1];
rz(-2.6182981) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6494203) q[0];
sx q[0];
rz(-1.5086552) q[0];
sx q[0];
rz(2.5323575) q[0];
rz(1.6388355) q[2];
sx q[2];
rz(-0.57300742) q[2];
sx q[2];
rz(-0.21493658) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1112422) q[1];
sx q[1];
rz(-2.6551464) q[1];
sx q[1];
rz(-0.72279795) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2449805) q[3];
sx q[3];
rz(-1.7953201) q[3];
sx q[3];
rz(2.6905439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3045197) q[2];
sx q[2];
rz(-2.4261116) q[2];
sx q[2];
rz(0.26930299) q[2];
rz(2.6473911) q[3];
sx q[3];
rz(-2.2952357) q[3];
sx q[3];
rz(2.0555029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.8158648) q[0];
sx q[0];
rz(-1.6115191) q[0];
sx q[0];
rz(1.4900526) q[0];
rz(1.4670463) q[1];
sx q[1];
rz(-2.84927) q[1];
sx q[1];
rz(-1.897859) q[1];
rz(-1.3080296) q[2];
sx q[2];
rz(-1.5599712) q[2];
sx q[2];
rz(0.47906265) q[2];
rz(0.79694637) q[3];
sx q[3];
rz(-1.4016101) q[3];
sx q[3];
rz(-2.3102643) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
