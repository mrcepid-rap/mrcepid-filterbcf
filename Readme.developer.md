# mrcepid-filterbcf Developer Readme

<!--
TODO: Please edit this Readme.developer.md file to include information
for developers or advanced users, for example:

* Information about app internals and implementation details
* How to report bugs or contribute to development
-->

## Testing

Note for development: for some unit tests the test data might be living on GCloud (if it's very large).
If you would like to gain access, please contact Eugene Gardner.
If you already have access to GCloud, run this command in the `/test_data/` directory:

```
gcloud storage rsync gs://iiuk-human-genetics/Syris/test_data/filterbcf .
```

## Running this app with additional computational resources

This app has the following entry points:

* main

When running this app, you can override the instance type to be used by
providing the ``systemRequirements`` field to ```/applet-XXXX/run``` or
```/app-XXXX/run```, as follows:

    {
      systemRequirements: {
        "main": {"instanceType": "mem2_hdd2_x2"}
      },
      [...]
    }

See <a
href="https://documentation.dnanexus.com/developer/api/running-analyses/io-and-run-specifications#run-specification">Run
Specification</a> in the API documentation for more information about the
available instance types.

