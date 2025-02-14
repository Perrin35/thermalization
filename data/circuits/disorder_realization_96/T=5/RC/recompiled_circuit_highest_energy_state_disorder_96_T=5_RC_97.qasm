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
rz(0.30492914) q[0];
sx q[0];
rz(-0.066216901) q[0];
sx q[0];
rz(0.15596341) q[0];
rz(1.1671542) q[1];
sx q[1];
rz(-2.4894297) q[1];
sx q[1];
rz(-0.48643938) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4733361) q[0];
sx q[0];
rz(-1.0205246) q[0];
sx q[0];
rz(-0.86007287) q[0];
x q[1];
rz(-0.42066853) q[2];
sx q[2];
rz(-1.0887869) q[2];
sx q[2];
rz(1.9276305) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.26775441) q[1];
sx q[1];
rz(-0.43711409) q[1];
sx q[1];
rz(0.74437352) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8016889) q[3];
sx q[3];
rz(-0.75711119) q[3];
sx q[3];
rz(-0.40386236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.277694) q[2];
sx q[2];
rz(-0.70964491) q[2];
sx q[2];
rz(-2.7126183) q[2];
rz(-0.046836827) q[3];
sx q[3];
rz(-0.39188477) q[3];
sx q[3];
rz(0.61280167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(1.2752537) q[0];
sx q[0];
rz(-0.99228042) q[0];
sx q[0];
rz(-0.66542768) q[0];
rz(-2.0003419) q[1];
sx q[1];
rz(-1.3803866) q[1];
sx q[1];
rz(-0.64249396) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2739853) q[0];
sx q[0];
rz(-1.6216177) q[0];
sx q[0];
rz(-2.4073868) q[0];
x q[1];
rz(1.8694473) q[2];
sx q[2];
rz(-0.89787102) q[2];
sx q[2];
rz(0.81097865) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4591532) q[1];
sx q[1];
rz(-2.2687468) q[1];
sx q[1];
rz(1.6070915) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8591381) q[3];
sx q[3];
rz(-2.5817462) q[3];
sx q[3];
rz(-2.0324872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.82681727) q[2];
sx q[2];
rz(-2.3544957) q[2];
sx q[2];
rz(-0.55244201) q[2];
rz(1.6636482) q[3];
sx q[3];
rz(-1.2976357) q[3];
sx q[3];
rz(2.4655931) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33573547) q[0];
sx q[0];
rz(-0.72000802) q[0];
sx q[0];
rz(2.3190401) q[0];
rz(2.3303253) q[1];
sx q[1];
rz(-2.8370116) q[1];
sx q[1];
rz(-1.6792345) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050348076) q[0];
sx q[0];
rz(-2.0322662) q[0];
sx q[0];
rz(1.1450923) q[0];
rz(-2.7285568) q[2];
sx q[2];
rz(-1.5025284) q[2];
sx q[2];
rz(0.78536805) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2671859) q[1];
sx q[1];
rz(-1.864434) q[1];
sx q[1];
rz(1.0974357) q[1];
rz(-1.261928) q[3];
sx q[3];
rz(-1.8654685) q[3];
sx q[3];
rz(2.1536649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.86639577) q[2];
sx q[2];
rz(-0.82922816) q[2];
sx q[2];
rz(0.59494507) q[2];
rz(-2.7368937) q[3];
sx q[3];
rz(-2.0345104) q[3];
sx q[3];
rz(-1.8428724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73612708) q[0];
sx q[0];
rz(-1.8647702) q[0];
sx q[0];
rz(-2.8986616) q[0];
rz(2.2566707) q[1];
sx q[1];
rz(-2.4156069) q[1];
sx q[1];
rz(0.0028217908) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045319917) q[0];
sx q[0];
rz(-0.82240538) q[0];
sx q[0];
rz(-2.1069645) q[0];
rz(-pi) q[1];
rz(-0.94945939) q[2];
sx q[2];
rz(-1.1295801) q[2];
sx q[2];
rz(-0.37829933) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0012008) q[1];
sx q[1];
rz(-2.2691548) q[1];
sx q[1];
rz(-0.091940885) q[1];
rz(-0.61496998) q[3];
sx q[3];
rz(-1.0677862) q[3];
sx q[3];
rz(-0.78395432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.60488492) q[2];
sx q[2];
rz(-1.9900091) q[2];
sx q[2];
rz(0.73337698) q[2];
rz(2.4865161) q[3];
sx q[3];
rz(-2.8337182) q[3];
sx q[3];
rz(-0.56183279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(0.68375278) q[0];
sx q[0];
rz(-0.34128749) q[0];
sx q[0];
rz(-2.999268) q[0];
rz(1.3617474) q[1];
sx q[1];
rz(-2.7889377) q[1];
sx q[1];
rz(0.49837643) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.179006) q[0];
sx q[0];
rz(-1.5091584) q[0];
sx q[0];
rz(0.85718244) q[0];
rz(1.519573) q[2];
sx q[2];
rz(-0.81256676) q[2];
sx q[2];
rz(2.0902772) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.315334) q[1];
sx q[1];
rz(-2.3439601) q[1];
sx q[1];
rz(0.41528292) q[1];
rz(-pi) q[2];
rz(2.8041995) q[3];
sx q[3];
rz(-2.3327347) q[3];
sx q[3];
rz(-2.0615426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.132823) q[2];
sx q[2];
rz(-1.885773) q[2];
sx q[2];
rz(-0.69457561) q[2];
rz(0.83235598) q[3];
sx q[3];
rz(-1.1135626) q[3];
sx q[3];
rz(2.8913403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.71448) q[0];
sx q[0];
rz(-0.50943333) q[0];
sx q[0];
rz(2.4795649) q[0];
rz(1.9561249) q[1];
sx q[1];
rz(-1.0292425) q[1];
sx q[1];
rz(1.1659291) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2625354) q[0];
sx q[0];
rz(-1.7723119) q[0];
sx q[0];
rz(-2.461947) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60816216) q[2];
sx q[2];
rz(-1.6972622) q[2];
sx q[2];
rz(1.6236931) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.59938972) q[1];
sx q[1];
rz(-1.521061) q[1];
sx q[1];
rz(-1.4982144) q[1];
x q[2];
rz(2.1816129) q[3];
sx q[3];
rz(-1.0418384) q[3];
sx q[3];
rz(2.4774266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7621496) q[2];
sx q[2];
rz(-1.2078441) q[2];
sx q[2];
rz(-0.70972788) q[2];
rz(2.659667) q[3];
sx q[3];
rz(-2.66633) q[3];
sx q[3];
rz(-0.019439241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4703281) q[0];
sx q[0];
rz(-0.98099357) q[0];
sx q[0];
rz(1.574466) q[0];
rz(1.6471242) q[1];
sx q[1];
rz(-2.6946805) q[1];
sx q[1];
rz(-2.2692197) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14022889) q[0];
sx q[0];
rz(-1.8793545) q[0];
sx q[0];
rz(-0.17805992) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.383423) q[2];
sx q[2];
rz(-2.2843666) q[2];
sx q[2];
rz(3.0564412) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1885729) q[1];
sx q[1];
rz(-1.5404048) q[1];
sx q[1];
rz(1.6769483) q[1];
rz(-1.5386343) q[3];
sx q[3];
rz(-0.65493203) q[3];
sx q[3];
rz(-0.82603588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.81277728) q[2];
sx q[2];
rz(-2.2810563) q[2];
sx q[2];
rz(-2.013773) q[2];
rz(0.40337107) q[3];
sx q[3];
rz(-1.6962467) q[3];
sx q[3];
rz(-0.71322125) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23680747) q[0];
sx q[0];
rz(-1.9984364) q[0];
sx q[0];
rz(-1.2971725) q[0];
rz(-0.5212658) q[1];
sx q[1];
rz(-1.8788985) q[1];
sx q[1];
rz(-3.0788132) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.990898) q[0];
sx q[0];
rz(-2.6932062) q[0];
sx q[0];
rz(-1.4175116) q[0];
x q[1];
rz(-0.36022236) q[2];
sx q[2];
rz(-1.7884322) q[2];
sx q[2];
rz(-1.9955903) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1181953) q[1];
sx q[1];
rz(-0.28701008) q[1];
sx q[1];
rz(1.0402518) q[1];
x q[2];
rz(2.2388458) q[3];
sx q[3];
rz(-2.3760736) q[3];
sx q[3];
rz(-2.3428041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.89821833) q[2];
sx q[2];
rz(-0.90460193) q[2];
sx q[2];
rz(0.96837366) q[2];
rz(-2.6280256) q[3];
sx q[3];
rz(-1.3526724) q[3];
sx q[3];
rz(-0.33094049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53720713) q[0];
sx q[0];
rz(-2.4810915) q[0];
sx q[0];
rz(-2.3576417) q[0];
rz(-1.2046658) q[1];
sx q[1];
rz(-0.37310633) q[1];
sx q[1];
rz(-1.7663667) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2581552) q[0];
sx q[0];
rz(-2.7984507) q[0];
sx q[0];
rz(0.68455066) q[0];
rz(1.7259459) q[2];
sx q[2];
rz(-1.5439234) q[2];
sx q[2];
rz(0.43153119) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.8177796) q[1];
sx q[1];
rz(-0.97517562) q[1];
sx q[1];
rz(2.2525674) q[1];
rz(-1.0635942) q[3];
sx q[3];
rz(-2.7873899) q[3];
sx q[3];
rz(2.2627047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2969926) q[2];
sx q[2];
rz(-2.832909) q[2];
sx q[2];
rz(-2.6958579) q[2];
rz(-2.5388057) q[3];
sx q[3];
rz(-1.2647537) q[3];
sx q[3];
rz(2.2980237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26850253) q[0];
sx q[0];
rz(-2.6188445) q[0];
sx q[0];
rz(-0.52421808) q[0];
rz(-1.2007319) q[1];
sx q[1];
rz(-1.8914765) q[1];
sx q[1];
rz(-1.1842309) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55726526) q[0];
sx q[0];
rz(-1.3818437) q[0];
sx q[0];
rz(-2.7897631) q[0];
rz(-pi) q[1];
rz(-1.2549868) q[2];
sx q[2];
rz(-1.5060052) q[2];
sx q[2];
rz(-3.0805317) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5752856) q[1];
sx q[1];
rz(-0.13338365) q[1];
sx q[1];
rz(-3.0493983) q[1];
rz(3.0940787) q[3];
sx q[3];
rz(-0.62944747) q[3];
sx q[3];
rz(0.87241064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.21620096) q[2];
sx q[2];
rz(-2.6915458) q[2];
sx q[2];
rz(1.0518543) q[2];
rz(2.9205186) q[3];
sx q[3];
rz(-1.3007921) q[3];
sx q[3];
rz(-0.93824798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5939519) q[0];
sx q[0];
rz(-1.5008391) q[0];
sx q[0];
rz(0.048901625) q[0];
rz(2.7652057) q[1];
sx q[1];
rz(-2.2793437) q[1];
sx q[1];
rz(1.9752165) q[1];
rz(1.6025887) q[2];
sx q[2];
rz(-2.1406662) q[2];
sx q[2];
rz(-1.6458423) q[2];
rz(-2.0509649) q[3];
sx q[3];
rz(-0.91560293) q[3];
sx q[3];
rz(-0.56165725) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
