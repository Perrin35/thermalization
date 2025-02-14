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
rz(1.4003657) q[0];
sx q[0];
rz(-0.23973149) q[0];
sx q[0];
rz(-2.8170407) q[0];
rz(-0.95207721) q[1];
sx q[1];
rz(-0.2258741) q[1];
sx q[1];
rz(1.3083375) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1150779) q[0];
sx q[0];
rz(-0.9008207) q[0];
sx q[0];
rz(1.7448533) q[0];
x q[1];
rz(-0.4842437) q[2];
sx q[2];
rz(-1.6855557) q[2];
sx q[2];
rz(1.5028624) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2497325) q[1];
sx q[1];
rz(-1.3123871) q[1];
sx q[1];
rz(2.051146) q[1];
rz(-1.7115373) q[3];
sx q[3];
rz(-2.2651197) q[3];
sx q[3];
rz(1.3545881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7294881) q[2];
sx q[2];
rz(-1.8821913) q[2];
sx q[2];
rz(2.7903902) q[2];
rz(-1.8515733) q[3];
sx q[3];
rz(-2.0100287) q[3];
sx q[3];
rz(-1.9180408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.7555162) q[0];
sx q[0];
rz(-1.2048683) q[0];
sx q[0];
rz(-2.2112041) q[0];
rz(-1.3049841) q[1];
sx q[1];
rz(-0.84140673) q[1];
sx q[1];
rz(-0.60943857) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064548858) q[0];
sx q[0];
rz(-1.6739558) q[0];
sx q[0];
rz(1.9984238) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0514293) q[2];
sx q[2];
rz(-1.8675641) q[2];
sx q[2];
rz(-2.0903843) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7623065) q[1];
sx q[1];
rz(-1.2067144) q[1];
sx q[1];
rz(-1.9696139) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8279161) q[3];
sx q[3];
rz(-2.5559396) q[3];
sx q[3];
rz(0.7440834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.096752) q[2];
sx q[2];
rz(-0.98016206) q[2];
sx q[2];
rz(0.074782221) q[2];
rz(-0.52037248) q[3];
sx q[3];
rz(-2.4986391) q[3];
sx q[3];
rz(2.5274932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.60270131) q[0];
sx q[0];
rz(-0.7170054) q[0];
sx q[0];
rz(0.30211788) q[0];
rz(-0.27613861) q[1];
sx q[1];
rz(-1.8243676) q[1];
sx q[1];
rz(-1.4322697) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0574422) q[0];
sx q[0];
rz(-2.1085116) q[0];
sx q[0];
rz(0.44769635) q[0];
rz(1.8319857) q[2];
sx q[2];
rz(-2.1495594) q[2];
sx q[2];
rz(-0.28489339) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.36784962) q[1];
sx q[1];
rz(-2.5605695) q[1];
sx q[1];
rz(1.2817205) q[1];
x q[2];
rz(1.4301705) q[3];
sx q[3];
rz(-2.6717477) q[3];
sx q[3];
rz(1.2805366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8351195) q[2];
sx q[2];
rz(-0.46572954) q[2];
sx q[2];
rz(-1.845537) q[2];
rz(0.45201388) q[3];
sx q[3];
rz(-2.0343503) q[3];
sx q[3];
rz(-0.38255102) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52962676) q[0];
sx q[0];
rz(-1.1671678) q[0];
sx q[0];
rz(-1.3324598) q[0];
rz(1.0491071) q[1];
sx q[1];
rz(-2.0499947) q[1];
sx q[1];
rz(-1.5528991) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72399607) q[0];
sx q[0];
rz(-1.231047) q[0];
sx q[0];
rz(-1.2509565) q[0];
rz(0.013601549) q[2];
sx q[2];
rz(-2.0715393) q[2];
sx q[2];
rz(0.73464962) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9178333) q[1];
sx q[1];
rz(-2.4661015) q[1];
sx q[1];
rz(1.2509327) q[1];
rz(-pi) q[2];
rz(-2.6060206) q[3];
sx q[3];
rz(-2.2102997) q[3];
sx q[3];
rz(0.65062338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66854746) q[2];
sx q[2];
rz(-1.0205525) q[2];
sx q[2];
rz(2.891053) q[2];
rz(-0.9969095) q[3];
sx q[3];
rz(-2.0018115) q[3];
sx q[3];
rz(-2.7610049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32960358) q[0];
sx q[0];
rz(-2.6530837) q[0];
sx q[0];
rz(-1.7937775) q[0];
rz(0.40428058) q[1];
sx q[1];
rz(-2.3826022) q[1];
sx q[1];
rz(1.1291198) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92496678) q[0];
sx q[0];
rz(-2.2499086) q[0];
sx q[0];
rz(0.51640262) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8724832) q[2];
sx q[2];
rz(-0.92293149) q[2];
sx q[2];
rz(-1.8108167) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.68539219) q[1];
sx q[1];
rz(-0.90185748) q[1];
sx q[1];
rz(0.13369707) q[1];
x q[2];
rz(0.32788289) q[3];
sx q[3];
rz(-2.8103069) q[3];
sx q[3];
rz(2.1123909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79444844) q[2];
sx q[2];
rz(-2.2114387) q[2];
sx q[2];
rz(1.244119) q[2];
rz(2.0813023) q[3];
sx q[3];
rz(-0.62842193) q[3];
sx q[3];
rz(-2.8384143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65619549) q[0];
sx q[0];
rz(-3.0363016) q[0];
sx q[0];
rz(-0.42718497) q[0];
rz(0.034491388) q[1];
sx q[1];
rz(-1.5692312) q[1];
sx q[1];
rz(3.131386) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9968352) q[0];
sx q[0];
rz(-3.138665) q[0];
sx q[0];
rz(0.63514955) q[0];
rz(-2.2030866) q[2];
sx q[2];
rz(-1.7006497) q[2];
sx q[2];
rz(-1.0009032) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31406826) q[1];
sx q[1];
rz(-1.6351103) q[1];
sx q[1];
rz(-0.8839279) q[1];
rz(-2.5764546) q[3];
sx q[3];
rz(-2.0083154) q[3];
sx q[3];
rz(-2.0034307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.61937845) q[2];
sx q[2];
rz(-0.41316119) q[2];
sx q[2];
rz(-2.348032) q[2];
rz(0.86269745) q[3];
sx q[3];
rz(-1.5746652) q[3];
sx q[3];
rz(0.70393744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4392387) q[0];
sx q[0];
rz(-0.84504253) q[0];
sx q[0];
rz(1.1335565) q[0];
rz(2.0639065) q[1];
sx q[1];
rz(-0.51529854) q[1];
sx q[1];
rz(-0.30212197) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0691119) q[0];
sx q[0];
rz(-1.9078507) q[0];
sx q[0];
rz(-1.2500136) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1094735) q[2];
sx q[2];
rz(-1.0240842) q[2];
sx q[2];
rz(-2.6448768) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6985642) q[1];
sx q[1];
rz(-0.57293597) q[1];
sx q[1];
rz(-0.40614508) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4601806) q[3];
sx q[3];
rz(-2.3294994) q[3];
sx q[3];
rz(-0.2821058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9083531) q[2];
sx q[2];
rz(-0.88963228) q[2];
sx q[2];
rz(-1.7944149) q[2];
rz(-1.2498648) q[3];
sx q[3];
rz(-1.1976539) q[3];
sx q[3];
rz(0.10716042) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5938479) q[0];
sx q[0];
rz(-0.43314728) q[0];
sx q[0];
rz(-0.10840848) q[0];
rz(-2.0938865) q[1];
sx q[1];
rz(-0.81755081) q[1];
sx q[1];
rz(1.0188867) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2813346) q[0];
sx q[0];
rz(-1.0509914) q[0];
sx q[0];
rz(-2.5297574) q[0];
x q[1];
rz(0.76129976) q[2];
sx q[2];
rz(-1.5660962) q[2];
sx q[2];
rz(-1.0710623) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68542504) q[1];
sx q[1];
rz(-1.9677094) q[1];
sx q[1];
rz(-2.6422068) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4014428) q[3];
sx q[3];
rz(-1.4664081) q[3];
sx q[3];
rz(-1.1889282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.63742796) q[2];
sx q[2];
rz(-2.1647858) q[2];
sx q[2];
rz(0.43761474) q[2];
rz(0.49312433) q[3];
sx q[3];
rz(-2.2212641) q[3];
sx q[3];
rz(-1.5776618) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2757932) q[0];
sx q[0];
rz(-1.9075305) q[0];
sx q[0];
rz(-0.27780521) q[0];
rz(1.5996784) q[1];
sx q[1];
rz(-0.96822396) q[1];
sx q[1];
rz(-2.7154198) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8086373) q[0];
sx q[0];
rz(-1.5530968) q[0];
sx q[0];
rz(-1.5874392) q[0];
rz(-0.8502281) q[2];
sx q[2];
rz(-0.93932187) q[2];
sx q[2];
rz(0.15635083) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8746169) q[1];
sx q[1];
rz(-2.8966581) q[1];
sx q[1];
rz(2.0744978) q[1];
x q[2];
rz(2.3087193) q[3];
sx q[3];
rz(-2.1578721) q[3];
sx q[3];
rz(-0.47490109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.49491945) q[2];
sx q[2];
rz(-0.16725954) q[2];
sx q[2];
rz(-2.0885928) q[2];
rz(0.35887512) q[3];
sx q[3];
rz(-1.5263298) q[3];
sx q[3];
rz(2.0487962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.7786355) q[0];
sx q[0];
rz(-0.21720049) q[0];
sx q[0];
rz(-2.2208075) q[0];
rz(2.6329363) q[1];
sx q[1];
rz(-0.76728907) q[1];
sx q[1];
rz(2.7210534) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9284436) q[0];
sx q[0];
rz(-0.9125114) q[0];
sx q[0];
rz(2.4930605) q[0];
x q[1];
rz(2.2013046) q[2];
sx q[2];
rz(-1.3740525) q[2];
sx q[2];
rz(1.3374752) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1433761) q[1];
sx q[1];
rz(-2.3399379) q[1];
sx q[1];
rz(0.90513913) q[1];
rz(-pi) q[2];
rz(-3.0402571) q[3];
sx q[3];
rz(-0.88235117) q[3];
sx q[3];
rz(0.42185129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4113808) q[2];
sx q[2];
rz(-2.3513942) q[2];
sx q[2];
rz(-2.8249557) q[2];
rz(2.5824879) q[3];
sx q[3];
rz(-0.63566256) q[3];
sx q[3];
rz(-0.81812286) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.485514) q[0];
sx q[0];
rz(-2.06388) q[0];
sx q[0];
rz(-2.0583454) q[0];
rz(-2.6651233) q[1];
sx q[1];
rz(-1.0404027) q[1];
sx q[1];
rz(-1.7146005) q[1];
rz(2.3219288) q[2];
sx q[2];
rz(-1.507187) q[2];
sx q[2];
rz(1.0095467) q[2];
rz(1.9540167) q[3];
sx q[3];
rz(-1.6578309) q[3];
sx q[3];
rz(2.1014392) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
