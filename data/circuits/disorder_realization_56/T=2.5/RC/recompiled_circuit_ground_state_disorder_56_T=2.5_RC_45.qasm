OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.174515) q[0];
sx q[0];
rz(1.4580026) q[0];
sx q[0];
rz(9.7249029) q[0];
rz(1.1974273) q[1];
sx q[1];
rz(-1.6150885) q[1];
sx q[1];
rz(1.6405029) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7952657) q[0];
sx q[0];
rz(-3.1259968) q[0];
sx q[0];
rz(-2.8898284) q[0];
rz(-pi) q[1];
rz(-1.9674703) q[2];
sx q[2];
rz(-0.80840092) q[2];
sx q[2];
rz(-1.3371181) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0847124) q[1];
sx q[1];
rz(-0.69082971) q[1];
sx q[1];
rz(0.34617306) q[1];
rz(-pi) q[2];
rz(0.0048801076) q[3];
sx q[3];
rz(-0.80880755) q[3];
sx q[3];
rz(-0.67833662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.40452051) q[2];
sx q[2];
rz(-1.9998735) q[2];
sx q[2];
rz(0.70297757) q[2];
rz(-2.5718555) q[3];
sx q[3];
rz(-2.2629786) q[3];
sx q[3];
rz(-1.5574633) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0074145929) q[0];
sx q[0];
rz(-2.5226722) q[0];
sx q[0];
rz(1.2331569) q[0];
rz(-0.59457072) q[1];
sx q[1];
rz(-2.1023127) q[1];
sx q[1];
rz(-2.0436683) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1933408) q[0];
sx q[0];
rz(-1.3701212) q[0];
sx q[0];
rz(-0.06788036) q[0];
rz(-pi) q[1];
rz(-2.4929201) q[2];
sx q[2];
rz(-0.59174109) q[2];
sx q[2];
rz(1.6206738) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8658376) q[1];
sx q[1];
rz(-0.081457327) q[1];
sx q[1];
rz(-1.2553196) q[1];
rz(3.0497562) q[3];
sx q[3];
rz(-0.79376924) q[3];
sx q[3];
rz(2.5240999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7132831) q[2];
sx q[2];
rz(-1.8417336) q[2];
sx q[2];
rz(-1.606288) q[2];
rz(-2.4226268) q[3];
sx q[3];
rz(-0.9459559) q[3];
sx q[3];
rz(0.21397056) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83800256) q[0];
sx q[0];
rz(-0.8135697) q[0];
sx q[0];
rz(-2.549262) q[0];
rz(-1.2447641) q[1];
sx q[1];
rz(-2.0015621) q[1];
sx q[1];
rz(-2.9959784) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8728565) q[0];
sx q[0];
rz(-1.4352192) q[0];
sx q[0];
rz(-0.4879293) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4136366) q[2];
sx q[2];
rz(-0.69456911) q[2];
sx q[2];
rz(0.63169152) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.023886746) q[1];
sx q[1];
rz(-2.5080006) q[1];
sx q[1];
rz(-2.0633477) q[1];
rz(-pi) q[2];
rz(2.0592732) q[3];
sx q[3];
rz(-1.8104189) q[3];
sx q[3];
rz(0.75877305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.738872) q[2];
sx q[2];
rz(-1.0307743) q[2];
sx q[2];
rz(-1.9090451) q[2];
rz(2.9018719) q[3];
sx q[3];
rz(-1.0545878) q[3];
sx q[3];
rz(-0.48503748) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0204912) q[0];
sx q[0];
rz(-0.70206577) q[0];
sx q[0];
rz(0.030601587) q[0];
rz(2.5864511) q[1];
sx q[1];
rz(-1.6629013) q[1];
sx q[1];
rz(-2.0487002) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26390938) q[0];
sx q[0];
rz(-2.110024) q[0];
sx q[0];
rz(0.79099057) q[0];
rz(-pi) q[1];
rz(-1.6299155) q[2];
sx q[2];
rz(-2.335304) q[2];
sx q[2];
rz(1.8102243) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2177201) q[1];
sx q[1];
rz(-1.7752547) q[1];
sx q[1];
rz(2.2108498) q[1];
rz(0.84597702) q[3];
sx q[3];
rz(-1.032077) q[3];
sx q[3];
rz(1.908345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0577724) q[2];
sx q[2];
rz(-1.3670992) q[2];
sx q[2];
rz(-0.27285451) q[2];
rz(1.7504292) q[3];
sx q[3];
rz(-0.83943668) q[3];
sx q[3];
rz(2.1896037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4680173) q[0];
sx q[0];
rz(-1.8033569) q[0];
sx q[0];
rz(-1.2982298) q[0];
rz(-0.6908373) q[1];
sx q[1];
rz(-1.1996484) q[1];
sx q[1];
rz(1.615049) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012155449) q[0];
sx q[0];
rz(-0.94506663) q[0];
sx q[0];
rz(2.246494) q[0];
rz(-pi) q[1];
rz(1.2882907) q[2];
sx q[2];
rz(-1.2307248) q[2];
sx q[2];
rz(-1.9486537) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4936798) q[1];
sx q[1];
rz(-2.4402752) q[1];
sx q[1];
rz(1.6108684) q[1];
x q[2];
rz(1.869321) q[3];
sx q[3];
rz(-1.9179452) q[3];
sx q[3];
rz(0.35856465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4096628) q[2];
sx q[2];
rz(-2.3356428) q[2];
sx q[2];
rz(2.9742677) q[2];
rz(-1.3903728) q[3];
sx q[3];
rz(-2.6087587) q[3];
sx q[3];
rz(0.65756857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9996027) q[0];
sx q[0];
rz(-2.7775192) q[0];
sx q[0];
rz(-1.2784736) q[0];
rz(-1.4356042) q[1];
sx q[1];
rz(-1.4878788) q[1];
sx q[1];
rz(-2.7659168) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4091051) q[0];
sx q[0];
rz(-1.8408436) q[0];
sx q[0];
rz(2.9069515) q[0];
rz(-pi) q[1];
rz(2.2907545) q[2];
sx q[2];
rz(-0.22432835) q[2];
sx q[2];
rz(-1.6825324) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3842114) q[1];
sx q[1];
rz(-0.39296341) q[1];
sx q[1];
rz(-2.1860366) q[1];
rz(-1.3270946) q[3];
sx q[3];
rz(-0.82984314) q[3];
sx q[3];
rz(2.0654701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0659236) q[2];
sx q[2];
rz(-2.0772987) q[2];
sx q[2];
rz(-0.86090487) q[2];
rz(1.9979477) q[3];
sx q[3];
rz(-2.836561) q[3];
sx q[3];
rz(2.8633269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025573108) q[0];
sx q[0];
rz(-1.8208068) q[0];
sx q[0];
rz(-0.75468165) q[0];
rz(-2.2017551) q[1];
sx q[1];
rz(-2.266423) q[1];
sx q[1];
rz(-1.8584724) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59903501) q[0];
sx q[0];
rz(-1.577566) q[0];
sx q[0];
rz(1.2640796) q[0];
x q[1];
rz(-0.91249429) q[2];
sx q[2];
rz(-0.55710885) q[2];
sx q[2];
rz(2.0832286) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.922328) q[1];
sx q[1];
rz(-1.1891051) q[1];
sx q[1];
rz(-1.8846518) q[1];
rz(0.51605172) q[3];
sx q[3];
rz(-0.95041227) q[3];
sx q[3];
rz(2.82016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5963886) q[2];
sx q[2];
rz(-2.0241006) q[2];
sx q[2];
rz(-2.0580573) q[2];
rz(-0.74294535) q[3];
sx q[3];
rz(-1.6094306) q[3];
sx q[3];
rz(1.7090181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98899406) q[0];
sx q[0];
rz(-0.78482634) q[0];
sx q[0];
rz(2.3866744) q[0];
rz(-2.6424291) q[1];
sx q[1];
rz(-2.1776336) q[1];
sx q[1];
rz(0.69854936) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23125739) q[0];
sx q[0];
rz(-2.4485579) q[0];
sx q[0];
rz(0.49654754) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6766308) q[2];
sx q[2];
rz(-1.4413646) q[2];
sx q[2];
rz(1.429582) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4515775) q[1];
sx q[1];
rz(-1.5507586) q[1];
sx q[1];
rz(-1.0385787) q[1];
rz(-pi) q[2];
rz(-0.0011092862) q[3];
sx q[3];
rz(-1.1517164) q[3];
sx q[3];
rz(2.1622879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.90427202) q[2];
sx q[2];
rz(-2.1492683) q[2];
sx q[2];
rz(0.80238706) q[2];
rz(2.5896416) q[3];
sx q[3];
rz(-0.80058432) q[3];
sx q[3];
rz(-0.62057453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.9629843) q[0];
sx q[0];
rz(-1.4405788) q[0];
sx q[0];
rz(-2.0340023) q[0];
rz(-1.7604609) q[1];
sx q[1];
rz(-1.3637204) q[1];
sx q[1];
rz(1.6345056) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0213623) q[0];
sx q[0];
rz(-1.7108546) q[0];
sx q[0];
rz(2.7441478) q[0];
rz(-pi) q[1];
rz(0.92169833) q[2];
sx q[2];
rz(-0.49921152) q[2];
sx q[2];
rz(-1.8550903) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9642051) q[1];
sx q[1];
rz(-1.9329762) q[1];
sx q[1];
rz(-1.7262579) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82689311) q[3];
sx q[3];
rz(-1.6312459) q[3];
sx q[3];
rz(-1.7884487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0570809) q[2];
sx q[2];
rz(-2.9456186) q[2];
sx q[2];
rz(1.8692807) q[2];
rz(1.48014) q[3];
sx q[3];
rz(-1.1353761) q[3];
sx q[3];
rz(-0.60012668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8965974) q[0];
sx q[0];
rz(-2.5801881) q[0];
sx q[0];
rz(-1.2835314) q[0];
rz(0.01783477) q[1];
sx q[1];
rz(-0.63648883) q[1];
sx q[1];
rz(3.0160115) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20480072) q[0];
sx q[0];
rz(-1.556634) q[0];
sx q[0];
rz(1.0165748) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4634976) q[2];
sx q[2];
rz(-1.1934278) q[2];
sx q[2];
rz(1.8630149) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.24841845) q[1];
sx q[1];
rz(-1.1496953) q[1];
sx q[1];
rz(1.1014654) q[1];
x q[2];
rz(1.3036895) q[3];
sx q[3];
rz(-1.3767164) q[3];
sx q[3];
rz(-1.1052856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77832001) q[2];
sx q[2];
rz(-1.2691701) q[2];
sx q[2];
rz(-0.8030836) q[2];
rz(0.17656365) q[3];
sx q[3];
rz(-1.8162138) q[3];
sx q[3];
rz(0.77809063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17978996) q[0];
sx q[0];
rz(-0.50146865) q[0];
sx q[0];
rz(1.5487221) q[0];
rz(1.8190307) q[1];
sx q[1];
rz(-1.7772728) q[1];
sx q[1];
rz(-1.5900236) q[1];
rz(1.5079458) q[2];
sx q[2];
rz(-1.812915) q[2];
sx q[2];
rz(-1.6356638) q[2];
rz(-0.33752058) q[3];
sx q[3];
rz(-1.6277085) q[3];
sx q[3];
rz(1.7276702) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
