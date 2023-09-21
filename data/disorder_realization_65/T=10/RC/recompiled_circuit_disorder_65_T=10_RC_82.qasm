OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31057519) q[0];
sx q[0];
rz(-1.0520881) q[0];
sx q[0];
rz(1.6488099) q[0];
rz(-2.2812023) q[1];
sx q[1];
rz(-0.69505039) q[1];
sx q[1];
rz(1.0746497) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.933691) q[0];
sx q[0];
rz(-1.5760734) q[0];
sx q[0];
rz(-3.123377) q[0];
rz(2.9825748) q[2];
sx q[2];
rz(-0.79610014) q[2];
sx q[2];
rz(3.0167992) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1638012) q[1];
sx q[1];
rz(-1.692354) q[1];
sx q[1];
rz(1.4896643) q[1];
rz(1.1172469) q[3];
sx q[3];
rz(-2.0587066) q[3];
sx q[3];
rz(-0.32353668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.82912123) q[2];
sx q[2];
rz(-2.1437936) q[2];
sx q[2];
rz(1.5134229) q[2];
rz(1.1734022) q[3];
sx q[3];
rz(-1.4459926) q[3];
sx q[3];
rz(2.9287958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5728977) q[0];
sx q[0];
rz(-0.64210367) q[0];
sx q[0];
rz(1.9182385) q[0];
rz(2.9808295) q[1];
sx q[1];
rz(-1.5784966) q[1];
sx q[1];
rz(1.5011903) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26186865) q[0];
sx q[0];
rz(-0.12517087) q[0];
sx q[0];
rz(-2.9414888) q[0];
rz(-pi) q[1];
rz(-2.1378527) q[2];
sx q[2];
rz(-2.6008285) q[2];
sx q[2];
rz(2.8218249) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7306819) q[1];
sx q[1];
rz(-2.4446179) q[1];
sx q[1];
rz(-2.8721786) q[1];
rz(-pi) q[2];
rz(0.0060175671) q[3];
sx q[3];
rz(-1.9618417) q[3];
sx q[3];
rz(2.1276377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.286065) q[2];
sx q[2];
rz(-2.3748886) q[2];
sx q[2];
rz(1.6274874) q[2];
rz(1.9747915) q[3];
sx q[3];
rz(-1.5224001) q[3];
sx q[3];
rz(0.91439009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1353772) q[0];
sx q[0];
rz(-1.051798) q[0];
sx q[0];
rz(-1.0789543) q[0];
rz(1.9619933) q[1];
sx q[1];
rz(-1.0410573) q[1];
sx q[1];
rz(-1.6378145) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30341879) q[0];
sx q[0];
rz(-1.6359328) q[0];
sx q[0];
rz(1.4831878) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90942803) q[2];
sx q[2];
rz(-0.69790188) q[2];
sx q[2];
rz(0.94240377) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7392002) q[1];
sx q[1];
rz(-1.270073) q[1];
sx q[1];
rz(-1.5261569) q[1];
rz(2.2272439) q[3];
sx q[3];
rz(-1.4667061) q[3];
sx q[3];
rz(-0.57746938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.43061313) q[2];
sx q[2];
rz(-2.6617699) q[2];
sx q[2];
rz(-0.7437931) q[2];
rz(-2.8841694) q[3];
sx q[3];
rz(-2.0580723) q[3];
sx q[3];
rz(-2.553489) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1733615) q[0];
sx q[0];
rz(-3.0259202) q[0];
sx q[0];
rz(-1.051735) q[0];
rz(0.60032088) q[1];
sx q[1];
rz(-2.5225263) q[1];
sx q[1];
rz(-1.1490885) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1968397) q[0];
sx q[0];
rz(-1.6819994) q[0];
sx q[0];
rz(2.8993594) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7096268) q[2];
sx q[2];
rz(-2.1199806) q[2];
sx q[2];
rz(0.84749046) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.035605343) q[1];
sx q[1];
rz(-0.78250256) q[1];
sx q[1];
rz(1.993202) q[1];
rz(1.9150919) q[3];
sx q[3];
rz(-2.2831884) q[3];
sx q[3];
rz(0.32574367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2465308) q[2];
sx q[2];
rz(-1.6677413) q[2];
sx q[2];
rz(2.4728298) q[2];
rz(-1.2459922) q[3];
sx q[3];
rz(-0.99530882) q[3];
sx q[3];
rz(-0.016383735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-0.65288654) q[0];
sx q[0];
rz(-1.1504494) q[0];
sx q[0];
rz(2.0892129) q[0];
rz(1.2106238) q[1];
sx q[1];
rz(-1.4167507) q[1];
sx q[1];
rz(-2.2573684) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0563755) q[0];
sx q[0];
rz(-1.7612805) q[0];
sx q[0];
rz(-2.9604244) q[0];
x q[1];
rz(-0.40712936) q[2];
sx q[2];
rz(-1.9756769) q[2];
sx q[2];
rz(-1.7083573) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4864038) q[1];
sx q[1];
rz(-0.15895325) q[1];
sx q[1];
rz(-2.7093191) q[1];
x q[2];
rz(1.0377117) q[3];
sx q[3];
rz(-1.248705) q[3];
sx q[3];
rz(-2.0562293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23652442) q[2];
sx q[2];
rz(-0.69513598) q[2];
sx q[2];
rz(-2.0489676) q[2];
rz(0.24400273) q[3];
sx q[3];
rz(-2.2677393) q[3];
sx q[3];
rz(2.3905335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7364175) q[0];
sx q[0];
rz(-0.002679499) q[0];
sx q[0];
rz(1.8238235) q[0];
rz(2.4353943) q[1];
sx q[1];
rz(-1.6786007) q[1];
sx q[1];
rz(0.96907369) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3450619) q[0];
sx q[0];
rz(-0.11632761) q[0];
sx q[0];
rz(-1.9274812) q[0];
rz(-pi) q[1];
rz(-1.6386119) q[2];
sx q[2];
rz(-0.8578476) q[2];
sx q[2];
rz(0.80578795) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.59086696) q[1];
sx q[1];
rz(-2.496965) q[1];
sx q[1];
rz(-2.4025737) q[1];
rz(8*pi/11) q[3];
sx q[3];
rz(-0.65398765) q[3];
sx q[3];
rz(-1.6884782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6270854) q[2];
sx q[2];
rz(-0.58834499) q[2];
sx q[2];
rz(0.43760854) q[2];
rz(-3.0833516) q[3];
sx q[3];
rz(-2.2831459) q[3];
sx q[3];
rz(-1.5313914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83201927) q[0];
sx q[0];
rz(-2.11634) q[0];
sx q[0];
rz(1.7818041) q[0];
rz(1.6925905) q[1];
sx q[1];
rz(-2.2429008) q[1];
sx q[1];
rz(1.4356027) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3754417) q[0];
sx q[0];
rz(-1.0784082) q[0];
sx q[0];
rz(0.96813162) q[0];
x q[1];
rz(0.82201634) q[2];
sx q[2];
rz(-3.1146345) q[2];
sx q[2];
rz(0.4617304) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3064733) q[1];
sx q[1];
rz(-1.8715579) q[1];
sx q[1];
rz(-1.9496099) q[1];
rz(-pi) q[2];
x q[2];
rz(0.42922677) q[3];
sx q[3];
rz(-1.1042522) q[3];
sx q[3];
rz(0.45688094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.74063611) q[2];
sx q[2];
rz(-1.2601968) q[2];
sx q[2];
rz(2.9675972) q[2];
rz(1.0081572) q[3];
sx q[3];
rz(-2.5139136) q[3];
sx q[3];
rz(-2.984916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.290264) q[0];
sx q[0];
rz(-2.2785318) q[0];
sx q[0];
rz(3.0086349) q[0];
rz(-2.6538387) q[1];
sx q[1];
rz(-0.98633352) q[1];
sx q[1];
rz(-1.3652323) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93850219) q[0];
sx q[0];
rz(-1.4346968) q[0];
sx q[0];
rz(-3.0654728) q[0];
rz(3.119602) q[2];
sx q[2];
rz(-1.4430265) q[2];
sx q[2];
rz(-0.86554229) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8147162) q[1];
sx q[1];
rz(-1.5183581) q[1];
sx q[1];
rz(-0.53175064) q[1];
x q[2];
rz(-0.41994628) q[3];
sx q[3];
rz(-1.4010324) q[3];
sx q[3];
rz(1.9672293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0662971) q[2];
sx q[2];
rz(-1.2434881) q[2];
sx q[2];
rz(0.42262849) q[2];
rz(-0.95831174) q[3];
sx q[3];
rz(-0.13923968) q[3];
sx q[3];
rz(-2.9039834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.002710297) q[0];
sx q[0];
rz(-0.30160987) q[0];
sx q[0];
rz(0.31059206) q[0];
rz(0.25767913) q[1];
sx q[1];
rz(-1.5035166) q[1];
sx q[1];
rz(0.92393595) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5101178) q[0];
sx q[0];
rz(-0.81561136) q[0];
sx q[0];
rz(-0.60879137) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6529979) q[2];
sx q[2];
rz(-1.5813507) q[2];
sx q[2];
rz(2.7985364) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.98665806) q[1];
sx q[1];
rz(-1.0549874) q[1];
sx q[1];
rz(0.45706473) q[1];
x q[2];
rz(2.4273848) q[3];
sx q[3];
rz(-1.076477) q[3];
sx q[3];
rz(-2.2480272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6442287) q[2];
sx q[2];
rz(-2.6022537) q[2];
sx q[2];
rz(1.3941992) q[2];
rz(-1.9723069) q[3];
sx q[3];
rz(-1.0401657) q[3];
sx q[3];
rz(2.6031301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7681463) q[0];
sx q[0];
rz(-0.099530749) q[0];
sx q[0];
rz(1.753153) q[0];
rz(-0.0069847981) q[1];
sx q[1];
rz(-0.81022898) q[1];
sx q[1];
rz(-3.1034234) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5725806) q[0];
sx q[0];
rz(-0.60927143) q[0];
sx q[0];
rz(1.6174181) q[0];
rz(1.1881234) q[2];
sx q[2];
rz(-1.6560935) q[2];
sx q[2];
rz(-1.6558937) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3341748) q[1];
sx q[1];
rz(-1.6125229) q[1];
sx q[1];
rz(1.9445022) q[1];
x q[2];
rz(-0.90922728) q[3];
sx q[3];
rz(-2.0709166) q[3];
sx q[3];
rz(-1.7928894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1756211) q[2];
sx q[2];
rz(-1.1493382) q[2];
sx q[2];
rz(1.9434628) q[2];
rz(-2.0576058) q[3];
sx q[3];
rz(-2.4626648) q[3];
sx q[3];
rz(-0.23588022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9772298) q[0];
sx q[0];
rz(-1.9259763) q[0];
sx q[0];
rz(1.5019793) q[0];
rz(-1.4981131) q[1];
sx q[1];
rz(-2.5479981) q[1];
sx q[1];
rz(2.610276) q[1];
rz(0.11001982) q[2];
sx q[2];
rz(-1.6286055) q[2];
sx q[2];
rz(1.5946228) q[2];
rz(-0.86919541) q[3];
sx q[3];
rz(-1.3251318) q[3];
sx q[3];
rz(-1.3983923) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];